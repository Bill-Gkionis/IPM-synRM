"""What the femm_queue tests share: processes that hold the slot, and stand-in tasks for the pool.

This is a module, not part of the test file, for the reason femm_test_tools.py gives: Windows starts a worker by
*spawning* a fresh Python, and that Python has to import the function it runs. pytest's own test modules cannot be
imported that way; a plain module in tests/ can, because conftest.py puts tests/ on sys.path.

A *stand-in* is a small made-up program that behaves like a solver in the one way a test looks at: it holds the slot for
a while, hangs for good, starts a child and says its pid, or reports its priority. No FEMM and no real solver is ever
started by the femm_queue tests, so they can run while the PC does other work.
"""
from __future__ import annotations

import multiprocessing
import os
import subprocess
import sys
import threading
import time
from pathlib import Path

import femm_test_tools as tools

import femm_queue

CTX = multiprocessing.get_context("spawn")           # the only way Windows starts a process; spelled out so it reads the same anywhere
STAND_IN = str(Path(__file__).with_name("stand_in_solver.py"))      # the stand-in solver programs (a script of its own)


# --------------------------------------------------------------------------------------- is a process still there?
def alive(pid: int) -> bool:
    """True if a process with this pid is running (the same Windows process list the femm_session tests use)."""
    return tools.process_info(pid) is not None


def wait_until_gone(pids, seconds: float = 15.0) -> list[int]:
    """Wait until none of these processes runs any more; return the ones that still do after `seconds` (an empty list is success)."""
    deadline = time.monotonic() + seconds
    while True:
        left = [pid for pid in pids if alive(pid)]
        if not left or time.monotonic() >= deadline:
            return left
        time.sleep(0.05)


# --------------------------------------------------------------------------------------- processes a test may end
# The rule: a test ends only a process it started, and only through a handle it holds. A bare pid is just a number, and Windows may
# give a dead process's number to a stranger later; a handle keeps its process's number from being given away while it is open, so
# a handle can only ever reach the process it was opened on. (Looking by pid, as alive() does, is harmless: it only reads.)
_SYNCHRONIZE = 0x00100000
_PROCESS_TERMINATE = 0x0001
_PROCESS_QUERY_LIMITED_INFORMATION = 0x1000
_WAIT_OBJECT_0, _WAIT_TIMEOUT = 0x0, 0x102
_kernel32_calls = None


def _kernel32():
    """kernel32 with the four calls Owned makes declared (a handle is 64 bits wide; ctypes assumes 32 unless told)."""
    global _kernel32_calls
    if _kernel32_calls is None:
        import ctypes
        from ctypes import wintypes
        k = ctypes.WinDLL("kernel32", use_last_error=True)
        k.OpenProcess.argtypes, k.OpenProcess.restype = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD], wintypes.HANDLE
        k.WaitForSingleObject.argtypes, k.WaitForSingleObject.restype = [wintypes.HANDLE, wintypes.DWORD], wintypes.DWORD
        k.TerminateProcess.argtypes, k.TerminateProcess.restype = [wintypes.HANDLE, wintypes.UINT], wintypes.BOOL
        k.CloseHandle.argtypes, k.CloseHandle.restype = [wintypes.HANDLE], wintypes.BOOL
        k.GetProcessTimes.argtypes = [wintypes.HANDLE] + [ctypes.POINTER(wintypes.FILETIME)] * 4
        k.GetProcessTimes.restype = wintypes.BOOL
        _kernel32_calls = k
    return _kernel32_calls


def _handle_creation_time(handle) -> int:
    """When the process behind this handle was created, as a Windows FILETIME (100 ns steps): with its pid, a process's identity."""
    import ctypes
    from ctypes import wintypes
    created, exited, kernel, user = (wintypes.FILETIME() for _ in range(4))
    if not _kernel32().GetProcessTimes(handle, ctypes.byref(created), ctypes.byref(exited), ctypes.byref(kernel), ctypes.byref(user)):
        raise ctypes.WinError(ctypes.get_last_error())
    return (created.dwHighDateTime << 32) | created.dwLowDateTime


class Owned:
    """A handle on a process that the test knows is its own, opened at the moment the process announced itself (so it was running).

    running() asks the handle, wait_gone() waits on it, end() terminates through it. Nothing here takes a pid again after the
    handle is open, so none of it can reach a process that merely got a dead one's number.

    `created` closes the last gap. A pid is a process's name only while it lives; if the stand-in had died before the handle was
    opened and Windows had given its number to a stranger in that moment, the handle would be the stranger's. The stand-in therefore
    announces its creation time too (pid and creation time together are what identify a process), and the handle is refused unless
    the time read through it is the announced one.
    """

    def __init__(self, pid: int, created: int | None = None):
        k = _kernel32()
        self.pid = pid
        self.handle = k.OpenProcess(_SYNCHRONIZE | _PROCESS_TERMINATE | _PROCESS_QUERY_LIMITED_INFORMATION, False, pid)
        if not self.handle:
            raise ProcessLookupError(f"process {pid} could not be opened (it had already ended)")
        if not self.running():
            self.close()
            raise ProcessLookupError(f"process {pid} had already ended when its handle was opened")
        self.created = _handle_creation_time(self.handle)
        if created is not None and self.created != created:
            self.close()
            raise ProcessLookupError(f"process {pid} is not the one that announced itself (its creation time differs: its number was reused)")

    def running(self) -> bool:
        return self.handle is not None and _kernel32().WaitForSingleObject(self.handle, 0) == _WAIT_TIMEOUT

    def wait_gone(self, seconds: float) -> bool:
        """True once the process has ended (within `seconds`)."""
        return self.handle is None or _kernel32().WaitForSingleObject(self.handle, int(seconds * 1000)) == _WAIT_OBJECT_0

    def end(self) -> None:
        if self.running():
            _kernel32().TerminateProcess(self.handle, 1)

    def close(self) -> None:
        if self.handle is not None:
            _kernel32().CloseHandle(self.handle)
            self.handle = None


def creation_time_of(pid: int) -> int:
    """The creation time of a running process. Only for a process that the caller holds (a pool's worker, say), so that the number is
    certainly its: the test that gets this value back can then check a handle against it (Owned(pid, created))."""
    owned = Owned(pid)
    try:
        return owned.created
    finally:
        owned.close()


class Tracker:
    """Opens an Owned handle on each stand-in process of a tree as it announces itself, and ends the leftovers through them.

    A stand-in announces itself by writing "<its pid> <its creation time>" into <folder>/<name>.pid (tests/stand_in_solver.py). The
    tracker polls for those files in a thread of its own, so a test that is blocked inside run_solver still has its handles in time,
    and it opens the handle at once, while the process is certainly running and certainly the test's own (it has just announced
    itself and is asleep), and refuses a handle whose creation time is not the announced one. The test then asks the handles whether
    anything survived, and on leaving the `with` the tracker ends whatever still runs, through those handles only, and closes them:
    a failed test leaves nothing running, and nothing else is ever touched.
    """

    def __init__(self, folder, names):
        self.folder, self.names = Path(folder), list(names)
        self.owned: dict[str, Owned | None] = {}               # None: announced, but gone before a handle could be opened
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._watch, daemon=True)

    def __enter__(self):
        self._thread.start()
        return self

    def _watch(self):
        while not self._stop.is_set() and len(self.owned) < len(self.names):
            for name in self.names:
                if name in self.owned:
                    continue
                try:
                    pid, created = (int(part) for part in (self.folder / f"{name}.pid").read_text().split())
                except (OSError, ValueError):                  # not there yet, or half written
                    continue
                try:
                    self.owned[name] = Owned(pid, created)
                except ProcessLookupError:
                    self.owned[name] = None
            time.sleep(0.02)

    def wait_announced(self, seconds: float = 60.0) -> bool:
        """True once every stand-in has announced itself and a handle is open on it."""
        deadline = time.monotonic() + seconds
        while len(self.owned) < len(self.names) and time.monotonic() < deadline:
            time.sleep(0.02)
        return len(self.owned) == len(self.names)

    def still_running(self, seconds: float = 15.0) -> list[str]:
        """The names of the processes that are still running `seconds` after this call (an empty list: all of them have ended)."""
        deadline = time.monotonic() + seconds
        left = []
        for name, owned in list(self.owned.items()):            # a copy: the watcher thread may still be adding to it
            if owned is not None and not owned.wait_gone(max(0.0, deadline - time.monotonic())):
                left.append(name)
        return left

    def __exit__(self, *exc):
        # A test that failed early may leave stand-ins that have not announced themselves yet (the program was started a moment ago):
        # give them a moment to do so, so that they too are held and ended. Bounded, and no delay at all when the test passed.
        self.wait_announced(5.0)
        self._stop.set()
        self._thread.join(10)
        for owned in self.owned.values():
            if owned is not None:
                owned.end()
                owned.close()


# --------------------------------------------------------------------------------------- a process that holds the slot
def _note(log, what: str) -> None:
    if log is not None:
        with open(log, "a") as handle:
            handle.write(f"{what} {time.time():.6f}\n")


def hold_slot(lock_dir, label, ready, release, log):
    """Runs in a process of its own: take the slot, say so, hold it until told to let go, then give it back.

    The "exit" note is written BEFORE the slot is given back, so the next holder's "enter" can never be earlier than it:
    that is what lets a test show that two holders never overlap.
    """
    with femm_queue.solver_slot(label, lock_dir=Path(lock_dir)):
        _note(log, "enter")
        ready.set()
        release.wait(120)
        _note(log, "exit")


class Holder:
    """Another process that holds the slot. Use it in a `with`: .pid is what holder.json must show, .release() lets go
    politely, .kill() ends it with no clean-up at all (what Task Manager's "End task" does)."""

    def __init__(self, lock_dir, label: str, log=None):
        self._ready, self._release = CTX.Event(), CTX.Event()
        self.process = CTX.Process(target=hold_slot, args=(str(lock_dir), label, self._ready, self._release,
                                                           None if log is None else str(log)))
        self.process.start()
        if not self._ready.wait(60):
            self.kill()
            raise AssertionError(f"the stand-in holder {label!r} did not get the slot within 60 s")

    @property
    def pid(self) -> int:
        return self.process.pid

    def release(self) -> None:
        self._release.set()
        self.process.join(30)
        assert self.process.exitcode == 0, f"the stand-in holder ended with exit code {self.process.exitcode}"

    def kill(self) -> None:
        if self.process.is_alive():
            self.process.kill()
        self.process.join(10)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        if self.process.is_alive():
            self._release.set()
            self.process.join(10)
        self.kill()


# --------------------------------------------------------------------------------------- a waiter that gets a Ctrl-C
def wait_then_get_interrupted(lock_dir, label, interrupt_after_s, outcome):
    """Runs in a process of its own: wait for a slot somebody else holds, and be interrupted after a while.

    signal.raise_signal(SIGINT) from a timer thread is what a Ctrl-C does to Python: the interpreter's own SIGINT handler
    runs, and it raises KeyboardInterrupt in the main thread, here inside femm_queue's waiting loop.
    """
    import signal
    import threading
    threading.Timer(interrupt_after_s, signal.raise_signal, args=(signal.SIGINT,)).start()
    try:
        with femm_queue.solver_slot(label, poll_s=0.2, lock_dir=Path(lock_dir)):
            outcome.put("got the slot")                       # must not happen: the holder never lets go
    except KeyboardInterrupt:
        outcome.put("interrupted")


def call_run_solver(lock_dir, folder, mode, seconds, timeout_s, interrupt_after_s, report):
    """Runs in a process of its own: run_solver on a stand-in program (tests/stand_in_solver.py, started in `mode`).

    With interrupt_after_s, a Ctrl-C arrives after that long, raised the way wait_then_get_interrupted raises it. What happened is put
    on `report`: ("returned", exit code, timed out), ("interrupted",), or ("raised", the error). A test that kills this process
    gets nothing on the report, which is the point of that test.
    """
    import signal
    if interrupt_after_s is not None:
        threading.Timer(interrupt_after_s, signal.raise_signal, args=(signal.SIGINT,)).start()
    try:
        result = femm_queue.run_solver([sys.executable, STAND_IN, mode, str(folder), str(seconds)], label="the caller", cwd=folder,
                                       timeout_s=timeout_s, lock_dir=Path(lock_dir))
        report.put(("returned", result.returncode, result.timed_out))
    except KeyboardInterrupt:
        report.put(("interrupted",))
    except BaseException as error:
        report.put(("raised", repr(error)))


# --------------------------------------------------------------------------------------- tasks for the pool
def echo_pid(item):
    """(the item, the pid of the worker that ran it)."""
    return item, os.getpid()


def maybe_stuck(args):
    """One task for the pool: (folder, index, mode). It writes one line per attempt, the pid of the worker that made it.

    mode "never": returns at once.  "first": the first attempt at this item hangs for good, a retry returns.
    "always": every attempt hangs. A hanging attempt first starts a child that also sleeps and writes its pid down:
    the stand-in for FEMM's fkn.exe, which lives below the worker and has to go with it.
    """
    folder, index, mode = args
    attempts = Path(folder) / f"item{index}.pids"
    with open(attempts, "a") as handle:
        handle.write(f"{os.getpid()}\n")
    attempt = len(attempts.read_text().split())
    if mode == "always" or (mode == "first" and attempt == 1):
        child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"])
        with open(Path(folder) / f"item{index}.children", "a") as handle:
            handle.write(f"{child.pid}\n")
        time.sleep(120)
    return index, os.getpid()


def hang_the_first_time(args):
    """One task for the pool: (folder, index, mode). It writes the pid of the worker that makes the attempt, one line a time, into <folder>/item<index>.pids.
    mode "first": the first attempt at this item hangs for good (a worker busy in a solve that never ends) and a retry returns; "never": it returns at once.
    Unlike maybe_stuck it starts no child of its own: what the worker started is its initializer's (a stand-in grandchild, see start_grandchild)."""
    folder, index, mode = args
    attempts = Path(folder) / f"item{index}.pids"
    with open(attempts, "a") as handle:
        handle.write(f"{os.getpid()}\n")
    if mode == "first" and len(attempts.read_text().split()) == 1:
        time.sleep(120)
    return index, os.getpid()


def maybe_die(args):
    """One task for the pool: (folder, index, mode), like maybe_stuck, but the attempt ends the whole worker with os._exit.

    That is a crash (a segfault in a library, the memory running out), not a hang. mode "first": only the first attempt crashes.
    "always": every attempt does.
    """
    folder, index, mode = args
    attempts = Path(folder) / f"item{index}.pids"
    with open(attempts, "a") as handle:
        handle.write(f"{os.getpid()}\n")
    attempt = len(attempts.read_text().split())
    if mode == "always" or (mode == "first" and attempt == 1):
        os._exit(3)
    return index, os.getpid()


def fail_on_three(item):
    if item == 3:
        raise ValueError("boom 3")
    return item


def missing_file(args):
    """A task that raises an OSError, the most ordinary error of a FEMM script (a model file that is not there). It writes one line
    per attempt, so a test can tell an error that was raised from a worker that was taken for dead and tried again."""
    folder, name = args
    with open(Path(folder) / "attempts.txt", "a") as handle:
        handle.write(f"{os.getpid()}\n")
    raise FileNotFoundError(f"no model file {name}")


def try_to_take_the_slot(lock_dir):
    """What a worker that tries to take the slot its own pool holds gets: a refusal (RuntimeError), not a wait."""
    try:
        with femm_queue.solver_slot("a worker", wait=False, lock_dir=Path(lock_dir)):
            return "TOOK THE SLOT"
    except femm_queue.QueueBusy as busy:
        return f"BUSY: {busy}"
    except RuntimeError as refused:
        return f"REFUSED: {refused}"


def hang_forever(pid_file):
    """An initializer that never returns: the stand-in for a FEMM that never answers while a worker starts."""
    Path(pid_file).write_text(str(os.getpid()))
    time.sleep(120)


def fail_to_start():
    raise OSError("no FEMM here (a stand-in for a FEMM that is not installed)")


def die_while_starting():
    """An initializer that ends the whole worker at once: a start-up crash, with no traceback to send."""
    os._exit(5)


def sleep_long(item):
    """A task that is still running two minutes from now (far longer than any test waits): a worker busy in a long solve."""
    time.sleep(120)


def pool_owner(lock_dir, pids_out):
    """Runs in a process of its own: own a pool of two workers, give each a task that does not end, report them.

    What is reported is a list of (pid, creation time) pairs, taken here, where the pool holds the workers (so each number is
    certainly that worker's), for the test to open handles on (Owned) and check against. The owner then sits inside map() until it
    is killed. The workers are BUSY, not idle, which is the case that matters: an idle worker notices a dead owner by itself (its
    pipe closes), a busy one is deep inside its task and would not.
    """
    with femm_queue.femm_pool("owner", workers=2, femm=False, lock_dir=Path(lock_dir)) as pool:
        pids_out.put([(pid, creation_time_of(pid)) for pid in pool.pids])
        pool.map(sleep_long, [0, 1])


# --------------------------------------------------------------------------------------- tasks for a pool with real FEMMs
def record_femm(folder):
    """An initializer. It runs after the worker's FEMM is open, so the binding already has a session: write down what it saw."""
    import femm
    import femm_session
    seen = femm._default is not None and femm_session.pid(femm._default) is not None
    Path(folder, f"initializer_{os.getpid()}.txt").write_text("femm open" if seen else "NO FEMM YET")


def solve_on_own_femm(args):
    """A task for a pool with FEMM: solve the learn/03 magnet on THIS worker's own FEMM. Returns (worker pid, FEMM pid, By in T).

    The half-second sleep keeps this worker busy long enough that the pool hands the next task to the other, idle, worker.
    """
    import femm
    import femm_session
    fem_path, hc = args
    time.sleep(0.5)
    by = tools.solve_magnet(fem_path, hc)
    return os.getpid(), femm_session.pid(femm._default), by



# --------------------------------------------------------------------------------------- the pool workers' jobs (the finding of a timeout trial of B9)
# A worker that is ended with no clean-up (TerminateProcess, End task, a crash, the pool ending it with Process.kill) must take everything it started with it. The tests start a
# stand-in grandchild below a worker and end the worker through a handle; the grandchild has to go with it. Nothing here starts a FEMM.
# The grandchild is the interpreter ITSELF, not sys.executable: in a virtual environment sys.executable is a small launcher in front of the real interpreter (two
# processes, two pids: the one a test would end is the launcher, and the interpreter would sleep on). multiprocessing starts its workers from sys._base_executable
# for the same reason. One process, one pid, one creation time: what a handle can be checked against.
PYTHON_ALONE = str(Path(sys.base_prefix) / "python.exe") if (Path(sys.base_prefix) / "python.exe").is_file() else sys.executable
SLEEP_TWO_MINUTES = [PYTHON_ALONE, "-I", "-c", "import time; time.sleep(120)"]      # the stand-in grandchild
SLEEP_A_MINUTE = [PYTHON_ALONE, "-I", "-c", "import time; time.sleep(60)"]          # the same for the worker that starts many: one that nobody records ends by itself soon


def announce(folder, name: str) -> None:
    """Write "<pid> <creation time>" of THIS process into <folder>/<name>_<pid>.txt: what a test needs to open a handle on exactly this process."""
    me = os.getpid()
    Path(folder, f"{name}_{me}.txt").write_text(f"{me} {creation_time_of(me)}")


def _record_grandchild(folder, child) -> None:
    with open(Path(folder, f"grandchildren_{os.getpid()}.txt"), "a") as handle:
        handle.write(f"{child.pid} {creation_time_of(child.pid)}\n")


def start_grandchild(folder, name: str = "worker"):
    """An initializer (it runs in every worker once, with no FEMM): announce this process, start the stand-in grandchild and write its
    "<pid> <creation time>" into <folder>/grandchildren_<this pid>.txt."""
    announce(folder, name)
    _record_grandchild(folder, subprocess.Popen(SLEEP_TWO_MINUTES))


def spawn_grandchildren(folder, count: int = 25, pause_s: float = 0.02):
    """An initializer that keeps this worker's MAIN thread starting processes for about a second: the situation that trial met, where a FEMM is
    created by the worker's main thread (open_femm) while the worker's owner is being ended. Every grandchild is written down as soon as it exists."""
    announce(folder, "worker")
    for _ in range(count):
        _record_grandchild(folder, subprocess.Popen(SLEEP_A_MINUTE))
        time.sleep(pause_s)


def own_job_facts(item):
    """A task for the pool: what this worker says about its own job (Windows). (its pid, femm_queue has made one, this process is in it,
    its handle is inheritable)."""
    import ctypes
    from ctypes import wintypes
    k = ctypes.WinDLL("kernel32", use_last_error=True)
    k.GetCurrentProcess.restype = wintypes.HANDLE
    k.IsProcessInJob.argtypes = [wintypes.HANDLE, wintypes.HANDLE, ctypes.POINTER(wintypes.BOOL)]
    k.IsProcessInJob.restype = wintypes.BOOL
    k.GetHandleInformation.argtypes = [wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD)]
    k.GetHandleInformation.restype = wintypes.BOOL
    job = getattr(femm_queue, "_OWN_JOB", None)
    if job is None:
        return os.getpid(), False, False, False
    inside, flags = wintypes.BOOL(), wintypes.DWORD()
    k.IsProcessInJob(k.GetCurrentProcess(), job.handle, ctypes.byref(inside))
    k.GetHandleInformation(job.handle, ctypes.byref(flags))
    return os.getpid(), True, bool(inside.value), bool(flags.value & 0x1)             # HANDLE_FLAG_INHERIT


def solve_and_count_the_job(args):
    """A task for a pool with FEMM: solve the learn/03 magnet on THIS worker's own FEMM, then say how many processes this worker's job has ever had
    (the worker itself, its FEMM, FEMM's fkn.exe of the solve) and how many it has now. Returns (worker pid, FEMM pid, By in T, total, active)."""
    import femm
    import femm_session
    fem_path, hc = args
    by = tools.solve_magnet(fem_path, hc)
    w = femm_queue._win()
    info = w.Accounting()
    if not w.k.QueryInformationJobObject(femm_queue._OWN_JOB.handle, 1, w.ctypes.byref(info), w.ctypes.sizeof(info), None):
        raise OSError("cannot read the accounting of this worker's job")
    return os.getpid(), femm_session.pid(femm._default), by, info.TotalProcesses, info.ActiveProcesses
