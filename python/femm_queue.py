"""femm_queue: one heavy solver run at a time on this PC, and a pool of FEMM workers that cannot hang it.

Why it exists
    A design rule (R10 of the tools plan): one heavy run at a time, machine-wide. Several sessions (a person at the desk, and the automated runs that work
    for that person) start solvers on the same PC: FEMM pools, CalculiX, Elmer, MATLAB, ROSS, SALib sweeps. Two heavy runs at once do
    not finish twice as fast; they slow each other down and the desktop stutters (ten FEMMs did it on 2026-10-28). The sessions
    cannot see each other's processes, so the rule has to live in code: every heavy run takes the ONE slot, runs at
    below-normal priority under a timeout, and can end nothing but the processes it started itself.
    Not covered yet (R8 of the tools plan): the study scripts that were there before this module (barrier_study.py, validate_full.py,
    optimize_demo.py, the flux_explorer sweeps, mesh_study.py, evaluate.py, worker_scaling.py, check_femm.py and the like) start their own
    FEMMs without the slot, and so does validation/run_all.py, which runs most of them in turn. Until the Electromagnetics session moves
    them onto femm_pool, run one of those only while nothing else heavy runs.

What is in it
    solver_slot(label)      take the slot (wait while somebody holds it) and keep it for the `with` block.
                            femm_slot is the first name of the same thing: one lock, not two
    holder()                who holds the slot now: {"label", "pid", "since"}, or None when it is free
    femm_pool(label, ...)   the slot, plus up to 8 worker processes with a FEMM each. pool.map(fn, items) -> list. A task over
                            its time has its worker's whole process tree ended and is tried once more
    run_solver(cmd, ...)    the slot, plus one program (ccx, ElmerSolver, ...) at below-normal priority under a timeout, in a Windows
                            job that the system ends with its caller
    QueueBusy               what wait=False raises when the slot is taken; its message names the holder
    python femm_queue.py    prints who holds the slot (exit status 0 free, 1 held) and takes nothing

How a script uses it
    import femm_queue

    # a FEMM batch: up to 8 workers, each with a FEMM of its own; a task over 120 s is ended and tried once more
    if __name__ == "__main__":                              # needed: every worker is a new Python that imports the script again
        with femm_queue.femm_pool("barrier study", task_timeout_s=120) as pool:
            rows = pool.map(solve_one_candidate, candidates)    # a module-level function, which a worker imports by name
    # a row that is {"status": "timeout"} is a candidate that hung twice; the rest of the batch is not affected

    # any other solver: one program, below normal priority, at most 4 threads, ended after 10 minutes
    run = femm_queue.run_solver([paths.tool("ccx"), "ring"], label="bench B4 ccx", cwd=job_folder, timeout_s=600, threads=4)
    run.returncode, run.stdout, run.stderr, run.seconds     # run.timed_out is True (and returncode None) if it was ended

    # only to look, or to fail at once instead of waiting
    femm_queue.holder()                                     # None, or {"label": ..., "pid": ..., "since": ...}
    with femm_queue.solver_slot("quick check", wait=False): ...    # raises QueueBusy, naming the holder, if it is taken

The lock, for a C programmer who is new to Python
    Two programs that must not run together need a mutex that both can see. A mutex lives in one process's memory, so
    they need something the operating system keeps for them. The OS keeps one for files: a process opens a file, then asks the OS to
    LOCK a few of its bytes. While that process holds the lock, nobody else's attempt to lock the same bytes succeeds. What the file
    contains does not matter; it is only the place where the OS keeps the lock, like a hook on a wall to hang a sign on.

    msvcrt.locking(fd, mode, nbytes) is Python's name for the C runtime's _locking() (Windows' LockFile underneath).
    `fd` is a file descriptor from os.open(), exactly as in C. The locked range starts at the CURRENT position of the file
    (that is why _try_lock() does os.lseek(fd, 0, ...) first, and so does _unlock()) and is nbytes long: one byte here.
    Mode LK_NBLCK means "lock, and if somebody else holds it, fail at once": Python raises PermissionError (errno EACCES).
    LK_UNLCK lets go. On Linux and macOS the same idea is fcntl.flock(fd, LOCK_EX | LOCK_NB), which the other branch of
    _try_lock() uses (that branch has not been run on this PC; the tests ran on Windows).

    Four facts about it carry the whole design (the first three were checked on this PC, 2026-11-05):
    1. The lock belongs to the OPEN FILE, the handle behind the fd, not to the file's name and not to the process. Two handles
       to the same file compete even inside one process. That is how holder() can ask "is it held?" by trying to take it.
    2. The OS lets go of a lock when the handle is closed, and it closes every handle of a process when the process ends, however
       it ends: finished, crashed, End task, taskkill /F. Nobody has to clean up. A holder killed with TerminateProcess left the
       lock free the next instant (tests/test_femm_queue.py: test_killed_holder_frees_the_queue).
    3. On Windows the lock is mandatory: while a byte is locked, another handle cannot even READ it (PermissionError). So
       the holder's name is not written into the lock file; it goes in holder.json beside it, where anybody can read it.
    4. The lock file is never deleted. If it were, a process that already has the old file open and one that creates a new file
       of the same name would each lock "the lock" and both would run.

    Why a lock and not a PID file? The C classic is "write my pid into x.pid, and if x.pid exists somebody is running". It breaks
    in two ways that a lock cannot. If the owner crashes, the file stays and says "running" for ever, and the next run has to
    GUESS whether pid N is still alive. And a number can come back: Windows is free to hand a dead process's number to a new,
    unrelated process later (not at once, but it does), so the file may then name a program that has nothing to do with the run and
    is alive. (Creating the file with O_EXCL makes the check-and-create atomic, so a race between two runs is solvable. Staleness
    and pid reuse are not.) With a lock "is it held?" is a question the kernel answers, atomically, about
    handles that are open right now; a process that is gone has none. holder.json is only a note for people, which is why
    holder() asks the lock first and reads the note second: a note left by a killed run is ignored.

Waiting, and what a waiter prints
    A caller that finds the slot taken prints one line per poll (every poll_s seconds, 30 by default):
        waiting for the solver slot: held by 'bench B4 ccx' (pid 12345) since 14:02
    It looks at the lock about once a second (not once per poll), so it starts within a second of the holder letting go. There is no
    queue order: whoever looks first after the release gets it, which is fair enough for runs of minutes. Ctrl-C while waiting
    ends the wait with KeyboardInterrupt and leaves nothing behind: nothing was taken, and the handle is closed.
    One block at a time per thread: a second solver_slot, femm_pool or run_solver in the thread that holds the slot would wait
    for itself for ever, so it raises RuntimeError instead. A worker of a pool, or a program run_solver started, inherits the
    environment variable MOTOR_SOLVER_SLOT ("<folder>|<holder pid>|<since>") and gets the same error. It is checked against
    the lock, so a stale copy of the variable (the holder is gone) is ignored. Pool tasks that run solvers themselves would be
    8 heavy runs, not one: start the solver from the pool's owner, not from a task.

The pool (femm_pool)
    Why processes, not threads: femm_session's docstring (one FEMM per process). Why not concurrent.futures.ProcessPoolExecutor:
    when one of its workers is killed it marks the WHOLE pool broken and every waiting task fails; multiprocessing.Pool replaces the
    worker but loses the task it was running, and map() never returns. So the pool here is small and its own. Each worker has a
    pipe to its owner. The owner's loop (Pool.map) hands a task to an idle worker, then waits for whichever comes first: a result,
    a worker that died, or a task's deadline. Nothing polls.
    A worker starts like this, in this order: lower its own priority; start a thread that ends the worker if the owner disappears, so a
    killed owner cannot leave workers computing after the slot is free (the thread starts no process); on Windows, join a JOB of its
    own (the kind run_solver's programs live in: kill every member when the last handle is closed, hold every member at the pool's
    priority class) before it starts any process, so that everything it starts afterwards, FEMM and its fkn.exe included, is born in
    it; femm_session.open_femm() (the priority is lowered before FEMM starts, so FEMM inherits it); then the caller's
    `initializer(*initargs)`. At its end it close_femm()s. The worker keeps the job's only handle, which is not inheritable, for
    its whole life (closing it would end the worker), so when the worker ends, however it ends (the pool ending it,
    TerminateProcess, End task, a crash, its own watch thread), the system closes the handle and ends every process the worker
    started. A worker that cannot join its job does not start FEMM or the caller's initializer: it tells its owner ("boot-failed", with the
    cause) and the pool refuses to open, the rule of run_solver: nothing is run unprotected, silently. (Why a job and not a walk of
    the worker's tree by the watch thread: a walk finds what exists when it looks, and a FEMM that the worker's main thread was just
    starting could be created after the walk and before the worker was gone. One was, in the timeout trial of B9 on 2026-11-06.)
    A task over task_timeout_s has its worker ended, and with it everything the worker started:
        Windows     the pool ends the worker through the handle it holds on it (Process.kill, which is TerminateProcess on that handle),
                    and the worker's own job ends the rest: FEMM 5.0 is the worker's child and starts a short-lived fkn.exe for each
                    solve, and all of them are members. No taskkill is started. taskkill /T walks the tree through each process's
                    RECORDED parent number, so an orphan whose dead parent's number the worker now holds could be taken for the worker's
                    child: a process this run did not start. A handle can only be the process it was opened on (Windows does not reuse
                    a number while a handle to the process is open), and a job holds exactly what the worker started.
        elsewhere   the worker leads its own process group, and the pool sends SIGKILL to the group (_end_tree).
    Never by name (a kill by image name ends every process of that name, an interactive FEMM window included). A fresh worker replaces the
    ended one at once, the task is tried once more, and a second timeout makes the item's result {"status": "timeout"} while the run
    goes on. A worker that dies by itself is handled the same way (the result is {"status": "died"}). An exception in a task is not a
    timeout: map() raises it, with the worker's traceback as its __cause__. A worker that does not start within start_timeout_s (a
    FEMM that never answers: callfemm has no timeout of its own) or whose initializer fails makes the pool refuse to open, with its
    traceback.

run_solver
    Starts the program at below-normal priority (on Windows creationflags=BELOW_NORMAL_PRIORITY_CLASS, so it is born at that
    class, and its job holds every process in it at that class: JOB_OBJECT_LIMIT_PRIORITY_CLASS, so a process that asks for a
    higher class, or starts a child with a higher one, stays where it is; MATLAB's graphics hosts did one of the two before the limit),
    with OMP_NUM_THREADS and CalculiX's own CCX_NPROC_EQUATION_SOLVER and CCX_NPROC_STIFFNESS set to `threads` in the
    CHILD's environment only. threads=None means the configured maximum, so a caller that forgets is still capped. It waits for
    the program for at most timeout_s (the time spent waiting for the slot does not count), ends it if that is over, and
    returns SolverRun(returncode, stdout, stderr, seconds). A program that ends by itself with a failure is a normal return
    (returncode 3, say); one that was ended on timeout has returncode None and .timed_out True.

    On Windows the program lives in a JOB OBJECT. A job is a group of processes that the operating system keeps track of as one:
    every process a member starts becomes a member, a member cannot leave, and the whole group can be ended in one call. This job
    has the limit "kill every member when the last handle to the job is closed", and only the caller holds a handle. What that gives:
    * A caller that is ended without clean-up (End task, a crash) cannot leave the program behind. The system closes the caller's
      handles: the lock is released AND the job is closed, which ends the program and everything it started, so the next holder
      does not start on top of a solver that nobody watches. (taskkill /T cannot give this: somebody alive has to run it, and it
      finds children through their parent, which an orphan no longer has.) The two releases are not ordered: the program is
      ended within a moment of the slot becoming free, not before it.
    * The program is created SUSPENDED (CREATE_SUSPENDED, together with the below-normal class), put into the job, and only then
      resumed, so no process of it can be born outside the job. Resuming takes a Toolhelp snapshot to find the program's one
      thread and calls ResumeThread. Every call is documented; NtResumeProcess would be a single call, but it is not documented.
    * run_solver returns only when the job has no process left, not when the program's own process exits: a launcher that hands its
      work to a child and exits at once (matlab.exe starting MATLAB.exe is the case to expect) must not free the slot while the
      heavy process runs on. The exit code stays the program's own, `seconds` runs until the whole tree has ended, and timeout_s
      covers the whole tree. (A helper that a solver leaves running for ever below it would therefore keep run_solver waiting
      until the timeout, and be ended with it.)
    * On timeout, on Ctrl-C and on an error the whole job is ended with TerminateJobObject, orphans included.
    * If the program cannot be put into the job (say the caller's own job forbids nesting), it is ended before it has run
      anything, and RuntimeError names the cause. A program is never run unprotected, silently.
    Not on Windows: the program leads its own process group, and a timeout or Ctrl-C ends the group with SIGKILL, but a caller that
    is itself ended without clean-up can leave the program running (the lock is released and nothing ends the program).
    Pool workers are covered on every system: on Windows each worker lives in a job of its own, which the system ends with the worker,
    however it ends; elsewhere each worker leads its own process group and watches its owner, and ends the group when the owner is gone.

The limits come from config.toml (paths.py), not from this file
    [femm] workers = 8: at most this many pool workers, and at most this many solver threads (R10).
    [femm] priority = 'below_normal': the class of every worker and solver. 'idle' is accepted too; 'normal' is refused (R10).
    [paths] out: the lock lives in out/.femm/ (solver.lock, holder.json).
    "Machine-wide" therefore means: every process that reads this config.toml. A second copy of the project somewhere else on the PC
    has an out/ of its own, so a slot of its own, and the two would not see each other. (There is one copy; keep it that way.)

What this module never does
    It never ends a process by name, a FEMM it did not start, or anything on port 48015 (nothing here calls femm.openfemm()).
    It ends exactly two kinds of process: a pool worker it started (on Windows through the handle the pool holds on it, and the worker's
    own job ends what the worker started; elsewhere through the worker's process group), and the one program run_solver started
    (with everything that program started: through its job on Windows, through its process group elsewhere). It starts no taskkill.

Run
    uv run python femm_queue.py                          who holds the slot now
    uv run pytest tests/test_femm_queue.py -q            the tests: no FEMM, no real solver; they use stand-in processes
    uv run python check_queue.py                         the live check: two real FEMM pools and a run_solver at the same time
Expect, with nothing running:  "the solver slot is free" and exit status 0. While a run holds it:
    "the solver slot is held by 'barrier study' (pid 12345) since 14:02" and exit status 1.
    The tests: all pass (on 2026-11-06 in about 40 s; the two that start FEMM are marked `femm`). The live check: see
    check_queue.py's docstring (on 2026-11-05 the three runs took turns in about 9 s and it ended in "PASS").

Self-check (the answers follow).
  1. A run holds the slot and is ended from Task Manager. How long until another run can take it, and what had to be cleaned up?
  2. Why does holder() try to take the lock instead of just reading holder.json?
  3. Why is the holder's label in a second file and not in the lock file?
  4. femm_pool starts eight workers, each calling open_femm(). Why does the pool take the slot once, and not each worker?
  5. A pool task calls run_solver(). What happens, and why is that better than letting it wait?
  6. On Windows a stuck worker is ended through the handle the pool holds on it, not with `taskkill /T /PID <pid>`, and never with
     `taskkill /IM femm_headlessw.exe`. Why?
  7. The caller of run_solver is ended from Task Manager while ccx runs. What happens to ccx on Windows, and why could taskkill /T
     not do the same?
Answers.
  1. An instant. The operating system closes the dead process's handles and the lock goes with them; nothing is cleaned up. The note
     (holder.json) stays behind and is ignored, because holder() asks the lock first.
  2. A note outlives its writer. After a crash it still names a run that is gone, and the process number in it may now belong to an
     unrelated program. The lock cannot be stale: it exists only while a live handle holds it.
  3. Windows locks are mandatory: nobody else can read a locked byte. The note has to be somewhere a waiting run can read.
  4. The slot is for the RUN, and eight workers are one run. The pool's owner already holds it. If each worker took it too, all eight
     would wait for the owner, and the owner waits for them: a deadlock. (The same reason a worker is refused if it asks.)
  5. It gets a RuntimeError at once. The pool's owner holds the slot and is waiting for the task, so the call could only wait for
     its own parent; and eight tasks running solvers would be eight heavy runs, not one.
  6. A handle is the process it was opened on and no other: Windows does not give a number away while a handle to the process is
     open, so the pool cannot hit a stranger that was given a dead process's number. taskkill /T finds what is below a process through
     the parent number each process recorded, so an orphan whose dead parent's number the worker now holds could be taken for the
     worker's child: a process this run did not start. The worker's own job holds exactly what the worker started (FEMM and its
     fkn.exe) and ends it with the worker. /IM ends every process of that name on the PC, including a FEMM that somebody has open by hand and that
     this run never started.
  7. The system ends it. ccx and everything it started are in a job that only the caller holds a handle to; when the caller dies its
     handles are closed, and closing the last handle of the job kills every member. taskkill /T has to be run by somebody, and a
     caller that was ended without clean-up never runs it; it also finds processes through their parent, so an orphan is out of its reach.
"""
from __future__ import annotations

import errno
import json
import multiprocessing
import multiprocessing.connection
import os
import pickle
import signal
import subprocess
import sys
import threading
import time
import traceback
from collections import deque
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import NamedTuple

import paths                                    # config.toml: [paths] out, [femm] workers and priority

if sys.platform == "win32":
    import msvcrt
else:
    import fcntl

__all__ = ["QueueBusy", "solver_slot", "femm_slot", "holder", "femm_pool", "Pool", "run_solver", "SolverRun",
           "WorkerTraceback", "main"]

LOCK_NAME = "solver.lock"                       # the file the operating system's lock sits on: never read, written or deleted
NOTE_NAME = "holder.json"                       # the note beside it that says who holds the lock
_CLAIM = "MOTOR_SOLVER_SLOT"                    # environment variable "<lock folder>|<holder pid>|<since>", inherited by children
_CTX = multiprocessing.get_context("spawn")     # how Windows starts a worker (a fresh Python); spelled out so it is the same anywhere
_PRIORITIES = ("below_normal", "idle")          # what config.toml's [femm] priority may say; 'normal' is refused on purpose (R10)
_START_TIMEOUT_S = 120.0                        # a worker that is not ready after this long did not start


class QueueBusy(RuntimeError):
    """The solver slot is taken and the caller said wait=False. The message names the holder; .holder is its note (or None)."""

    def __init__(self, message: str, holder: dict | None = None):
        super().__init__(message)
        self.holder = holder


# ============================================================================================== the lock (the operating system's)
def _try_lock(fd: int) -> bool:
    """Try to take the OS lock on the file behind `fd`, without waiting. True: we hold it now. False: somebody else does."""
    if sys.platform == "win32":
        os.lseek(fd, 0, os.SEEK_SET)                    # msvcrt.locking locks the bytes from the CURRENT position
        try:
            msvcrt.locking(fd, msvcrt.LK_NBLCK, 1)      # one byte, exclusive, and do not wait
        except OSError as error:
            if error.errno in (errno.EACCES, errno.EDEADLK):    # "locking violation": another handle holds that byte
                return False
            raise                                       # anything else (bad fd, ...) is a bug, not a busy slot
        return True
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:                             # EWOULDBLOCK: another open file description holds it
        return False
    return True


def _unlock(fd: int) -> None:
    """Let go of the lock on `fd`. Closing the handle would do it too; doing it first makes the release immediate."""
    if sys.platform == "win32":
        os.lseek(fd, 0, os.SEEK_SET)                    # LK_UNLCK also works from the current position
        msvcrt.locking(fd, msvcrt.LK_UNLCK, 1)
    else:
        fcntl.flock(fd, fcntl.LOCK_UN)


# ============================================================================================== the note: who holds it
def _lock_folder(lock_dir) -> Path:
    """The folder of the lock: `lock_dir` if given (the tests do), else out/.femm of config.toml's [paths] out."""
    return Path(lock_dir) if lock_dir is not None else paths.path("out") / ".femm"


def _key(folder: Path) -> str:
    """A folder as one comparable string (Windows ignores case in paths)."""
    return os.path.normcase(os.path.abspath(folder))


def _hhmm(since) -> str:
    try:
        return datetime.fromisoformat(since).strftime("%H:%M")
    except (TypeError, ValueError):
        return "?"


def _who(info: dict | None) -> str:
    """"'label' (pid N) since HH:MM": how a holder is named in a waiting line and in QueueBusy's message."""
    if not info or info.get("pid") is None:
        return "a run whose note (holder.json) cannot be read"
    return f"'{info['label']}' (pid {info['pid']}) since {_hhmm(info['since'])}"


def _write_note(folder: Path, label: str) -> tuple[dict, bool]:
    """Write holder.json for this run. (the note, whether it could be written). The lock is what holds the slot; the note is
    only for people, so a failure to write it is a warning, not a reason to give the slot up."""
    info = {"label": label, "pid": os.getpid(), "since": datetime.now().astimezone().isoformat(timespec="seconds")}
    try:
        (folder / NOTE_NAME).write_text(json.dumps(info), encoding="utf-8")
        return info, True
    except OSError as error:
        print(f"warning: cannot write {folder / NOTE_NAME} ({error}): the slot is held, but a waiting run cannot name the holder", flush=True)
        try:
            (folder / NOTE_NAME).unlink()               # do not leave the note of a run that has ended in its place
        except OSError:
            pass
        return info, False


def _remove_note(folder: Path) -> None:
    try:
        (folder / NOTE_NAME).unlink()
    except OSError:                                     # already gone, or a reader has it open this instant: the note is only a note
        pass


def _read_note(folder: Path, patience_s: float = 0.5) -> dict | None:
    """holder.json as a dict, or None. Only called when the lock IS held, so a missing or half-written note means its
    writer is between taking the lock and writing the note: look again for up to `patience_s` before giving up."""
    deadline = time.monotonic() + patience_s
    while True:
        try:
            info = json.loads((folder / NOTE_NAME).read_text(encoding="utf-8"))
            if isinstance(info, dict) and {"label", "pid", "since"} <= set(info):
                return info
        except (OSError, ValueError):
            pass
        if time.monotonic() >= deadline:
            return None
        time.sleep(0.02)


def holder(lock_dir=None) -> dict | None:
    """Who holds the solver slot now: {"label", "pid", "since"} ("since" is an ISO time of day), or None if it is free.

    The LOCK decides, not the note: this tries to take the lock for an instant. If that works, nobody holds the slot, whatever
    holder.json says (a killed run leaves its note behind). If it does not work, somebody does, and the note says who. A
    holder that could not be named (its note is missing or unreadable) is reported as {"label": "unknown", "pid": None,
    "since": None}. Asking creates nothing.
    """
    folder = _lock_folder(lock_dir)
    lock_file = folder / LOCK_NAME
    if not lock_file.exists():
        return None
    fd = os.open(lock_file, os.O_RDWR)
    try:
        if _try_lock(fd):
            _unlock(fd)
            return None
    finally:
        os.close(fd)
    return _read_note(folder) or {"label": "unknown", "pid": None, "since": None}


# ============================================================================================== the slot
_local = threading.local()                      # per thread: the lock folders this thread holds the slot for


def _held_here() -> set:
    held = getattr(_local, "folders", None)
    if held is None:
        held = _local.folders = set()
    return held


def _refuse_nesting(folder: Path) -> None:
    """Raise RuntimeError, instead of waiting for ever, when this thread (or the run that started this process) holds the slot."""
    key = _key(folder)
    if key in _held_here():
        raise RuntimeError(
            "this thread already holds the solver slot (a solver_slot, femm_pool or run_solver block is still open here), so "
            "asking for it again would wait for itself for ever. Leave that block first, or take the slot once around the whole "
            "run: femm_pool and run_solver take it themselves, so do not wrap them in solver_slot.")
    claim = os.environ.get(_CLAIM)
    if not claim:
        return
    try:
        claimed_key, claimed_pid, claimed_since = claim.rsplit("|", 2)
        claimed_pid = int(claimed_pid)
    except ValueError:
        return
    if claimed_key != key or claimed_pid == os.getpid():        # another lock folder, or this very process (another thread)
        return
    info = holder(folder)                                       # the lock is the truth: is the claim still good?
    if info and info["pid"] == claimed_pid and info["since"] == claimed_since:
        raise RuntimeError(
            f"this process was started by a run that holds the solver slot ({_who(info)}), so asking for it here would wait for "
            "that run, which is waiting for this process. A worker or a solver started under the slot already runs under it. "
            "(If this is a pool worker that re-ran the script: a script that opens a pool needs the usual "
            "`if __name__ == \"__main__\":` guard, because every worker imports it again.)")


def _take(fd: int, folder: Path, wait: bool, poll_s: float) -> None:
    """Return once `fd` holds the lock. If it is taken: QueueBusy when wait is False, else print a line per poll and keep looking."""
    last_line = None
    while True:
        for _ in range(4):                                      # a few quick tries before calling it busy: holder() probes the lock
            if _try_lock(fd):                                   # for an instant, and a holder letting go takes about as long
                return
            time.sleep(0.01)
        info = _read_note(folder)
        if not wait:
            raise QueueBusy(f"the solver slot is held by {_who(info)}", info)
        now = time.monotonic()
        if last_line is None or now - last_line >= poll_s:
            print(f"waiting for the solver slot: held by {_who(info)}", flush=True)
            last_line = now
        time.sleep(min(poll_s, 1.0))                            # look at the lock about once a second, print once per poll


@contextmanager
def solver_slot(label: str, *, wait: bool = True, poll_s: float = 30.0, lock_dir=None):
    """Hold the machine-wide solver slot for the `with` block; yields the holder note {"label", "pid", "since"}.

    wait=True (the default): while somebody else holds it, print `waiting for the solver slot: held by '<label>' (pid N) since
    HH:MM` every poll_s seconds and keep looking. wait=False: raise QueueBusy, naming the holder, at once. The slot is given
    back however the block ends (return, exception, Ctrl-C), and by the operating system if the process is killed.
    RuntimeError, not a wait, if this thread already holds it or this process was started by the run that does (module docstring).
    lock_dir: the lock folder, by default out/.femm; the tests pass their own so that they never touch the real slot.
    """
    if not poll_s > 0:
        raise ValueError(f"poll_s={poll_s!r}: a positive number of seconds")
    folder = _lock_folder(lock_dir)
    _refuse_nesting(folder)
    folder.mkdir(parents=True, exist_ok=True)
    fd = os.open(folder / LOCK_NAME, os.O_RDWR | os.O_CREAT)    # not inheritable: no child process ever holds this handle
    locked = noted = False
    try:
        _take(fd, folder, wait, poll_s)
        locked = True
        info, noted = _write_note(folder, label)
        key, previous = _key(folder), os.environ.get(_CLAIM)
        try:                                                    # the markers are set INSIDE the try: an interrupt cannot leave them behind
            _held_here().add(key)
            os.environ[_CLAIM] = f"{key}|{info['pid']}|{info['since']}"
            yield info
        finally:
            _held_here().discard(key)
            if previous is None:
                os.environ.pop(_CLAIM, None)
            else:
                os.environ[_CLAIM] = previous
    finally:
        if noted:
            _remove_note(folder)                                # before the unlock: after it the note would be the next holder's
        if locked:
            try:
                _unlock(fd)
            except OSError:                                     # closing the handle below releases it anyway
                pass
        os.close(fd)


femm_slot = solver_slot                         # the first name (FEMM only); the same context manager, the same lock (R10)


# ============================================================================================== priority, and ending a process tree
def _priority_name() -> str:
    """config.toml's [femm] priority, refusing anything that is not below normal (R10)."""
    name = paths.setting("femm", "priority")
    if name not in _PRIORITIES:
        raise ValueError(f"[femm] priority in {paths.CONFIG} is {name!r}: solvers run at one of {', '.join(map(repr, _PRIORITIES))} "
                         "(rule R10: never at normal priority)")
    return name


def _limit(what: str, value, minimum: int = 1) -> int:
    """`value` as a count between `minimum` and the configured maximum ([femm] workers); None means the maximum."""
    cap = paths.setting("femm", "workers")
    if value is None:
        return cap
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= cap:
        raise ValueError(f"{what}={value!r}: a whole number between {minimum} and {cap} (config.toml [femm] workers, rule R10)")
    return value


def _windows_class(name: str) -> int:
    return {"below_normal": subprocess.BELOW_NORMAL_PRIORITY_CLASS, "idle": subprocess.IDLE_PRIORITY_CLASS}[name]


def _lower_self(name: str) -> None:
    """Put this process at the configured priority, never raising it (a process that is already idle stays idle)."""
    if sys.platform == "win32":
        import ctypes
        from ctypes import wintypes
        rank = {subprocess.IDLE_PRIORITY_CLASS: 0, subprocess.BELOW_NORMAL_PRIORITY_CLASS: 1, subprocess.NORMAL_PRIORITY_CLASS: 2,
                subprocess.ABOVE_NORMAL_PRIORITY_CLASS: 3, subprocess.HIGH_PRIORITY_CLASS: 4, 0x100: 5}
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetCurrentProcess.restype = wintypes.HANDLE
        kernel32.GetPriorityClass.argtypes = [wintypes.HANDLE]
        kernel32.GetPriorityClass.restype = wintypes.DWORD
        kernel32.SetPriorityClass.argtypes = [wintypes.HANDLE, wintypes.DWORD]
        kernel32.SetPriorityClass.restype = wintypes.BOOL
        me, wanted = kernel32.GetCurrentProcess(), _windows_class(name)
        if rank.get(kernel32.GetPriorityClass(me), 2) > rank[wanted]:
            if not kernel32.SetPriorityClass(me, wanted):
                raise ctypes.WinError(ctypes.get_last_error())
    else:
        target = {"below_normal": 5, "idle": 19}[name]
        niceness = os.getpriority(os.PRIO_PROCESS, 0)
        if niceness < target:
            os.nice(target - niceness)


def _end_tree(pid: int) -> tuple[int, str]:
    """Not on Windows: end the process group that `pid` leads, with SIGKILL. (0, "") if the signal was sent, (128, "no such process")
    or (1, the error's text) if not. BY PID, NEVER BY NAME.

    Only ever called with the pid of a process this module started, which leads its own session (start_new_session / setsid), so
    its process group is its pid: the program of run_solver (_end_program), a pool worker (Pool._end_worker) and, from a worker's own
    watch thread, the worker itself (_watch_owner).

    Windows: nothing calls it, and it has no taskkill branch any more. run_solver ends its program through its job; the pool ends a worker
    through the handle it holds on it (Process.kill) and the worker's own job ends what the worker started; the watch thread only ends its
    own worker (os._exit) and leaves the rest to that job. The first version walked the tree with taskkill /T /F /PID, which finds children
    through the parent number each process recorded: an orphan whose dead parent's number the worker now holds could be taken for the
    worker's child, a process this run did not start (the final review of the tools plan, Important 2). Called on Windows it raises.
    """
    if sys.platform == "win32":
        raise RuntimeError("_end_tree ends a process group and Windows has none: a pool worker is ended through its handle and its job, "
                           "and run_solver's program through its job")
    try:
        os.killpg(pid, signal.SIGKILL)
        return 0, ""
    except ProcessLookupError:
        return 128, "no such process"
    except PermissionError as error:
        return 1, str(error)


# ============================================================================================== run_solver
class SolverRun(NamedTuple):
    """What run_solver returns: the exit code (None if the program was ended on timeout), its output, and how long it ran."""
    returncode: int | None
    stdout: str
    stderr: str
    seconds: float

    @property
    def timed_out(self) -> bool:
        return self.returncode is None


def _popen_options(priority: str) -> dict:
    """Off Windows: how to start a program so that it is born at the configured priority and leads its own process group."""
    return {"start_new_session": True, "preexec_fn": lambda: os.nice({"below_normal": 5, "idle": 19}[priority])}


def _text(partial) -> str:
    if partial is None:
        return ""
    return partial.decode("utf-8", "replace") if isinstance(partial, bytes) else partial


def _communicate(proc: subprocess.Popen, deadline: float, timeout_s: float) -> tuple[str, str]:
    """proc.communicate() until `deadline` (a time.monotonic() value), in slices of at most half a second.

    One long communicate(timeout=...) cannot be interrupted by Ctrl-C on Windows: it waits in Thread.join(timeout), which a
    signal does not wake (measured on this PC, 2026-11-05: a Ctrl-C raised at 1 s was noticed at 30 s, when the wait ended).
    In slices a Ctrl-C is noticed within half a second. Raises subprocess.TimeoutExpired when the deadline passes.
    """
    while True:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise subprocess.TimeoutExpired(proc.args, timeout_s)
        try:
            return proc.communicate(timeout=min(remaining, 0.5))
        except subprocess.TimeoutExpired:
            continue


def _last_output(proc: subprocess.Popen) -> tuple[str, str]:
    """What the program wrote before it was ended (its pipes close once the whole tree is gone)."""
    try:
        return proc.communicate(timeout=10)
    except subprocess.TimeoutExpired as late:                   # something that escaped still holds the pipes: keep what was read
        return _text(late.stdout), _text(late.stderr)


def _reap(proc: subprocess.Popen) -> None:
    """After the program was ended on an interrupt or an error: collect its pipes and exit code. Never raises (but for Ctrl-C)."""
    try:
        proc.communicate(timeout=10)
    except Exception:
        pass


def _end_program(proc: subprocess.Popen) -> None:
    """Off Windows: end the program that run_solver started, with everything below it, and wait for it to be gone."""
    for _ in range(2):
        if proc.poll() is not None:
            return
        code, text = _end_tree(proc.pid)
        try:
            proc.wait(10)
            return
        except subprocess.TimeoutExpired:
            continue
    raise RuntimeError(f"could not end the program (pid {proc.pid}); it said {code}: {text}")


# -------------------------------------------------------------------------------------- Windows: the program's job object
_CREATE_SUSPENDED = 0x00000004                  # CreateProcess flag: the program's thread is created but does not run yet
_WIN = None                                     # kernel32 with the job calls declared, and the structures; made once by _win()


def _win():
    """kernel32 with the calls of the job code declared (argument and result types, as in femm_session._version_calls: without
    them ctypes treats every handle as a 32-bit int), and the structures they fill. Built once, on first use."""
    global _WIN
    if _WIN is None:
        import ctypes
        import types
        from ctypes import wintypes
        k = ctypes.WinDLL("kernel32", use_last_error=True)

        class BasicLimit(ctypes.Structure):             # JOBOBJECT_BASIC_LIMIT_INFORMATION
            _fields_ = [("PerProcessUserTimeLimit", ctypes.c_int64), ("PerJobUserTimeLimit", ctypes.c_int64),
                        ("LimitFlags", wintypes.DWORD), ("MinimumWorkingSetSize", ctypes.c_size_t),
                        ("MaximumWorkingSetSize", ctypes.c_size_t), ("ActiveProcessLimit", wintypes.DWORD),
                        ("Affinity", ctypes.c_size_t), ("PriorityClass", wintypes.DWORD), ("SchedulingClass", wintypes.DWORD)]

        class IoCounters(ctypes.Structure):             # IO_COUNTERS
            _fields_ = [(name, ctypes.c_uint64) for name in (
                "ReadOperationCount", "WriteOperationCount", "OtherOperationCount",
                "ReadTransferCount", "WriteTransferCount", "OtherTransferCount")]

        class ExtendedLimit(ctypes.Structure):          # JOBOBJECT_EXTENDED_LIMIT_INFORMATION (144 bytes on 64-bit Windows)
            _fields_ = [("BasicLimitInformation", BasicLimit), ("IoInfo", IoCounters), ("ProcessMemoryLimit", ctypes.c_size_t),
                        ("JobMemoryLimit", ctypes.c_size_t), ("PeakProcessMemoryUsed", ctypes.c_size_t),
                        ("PeakJobMemoryUsed", ctypes.c_size_t)]

        class Accounting(ctypes.Structure):             # JOBOBJECT_BASIC_ACCOUNTING_INFORMATION (48 bytes)
            _fields_ = [("TotalUserTime", ctypes.c_int64), ("TotalKernelTime", ctypes.c_int64),
                        ("ThisPeriodTotalUserTime", ctypes.c_int64), ("ThisPeriodTotalKernelTime", ctypes.c_int64),
                        ("TotalPageFaultCount", wintypes.DWORD), ("TotalProcesses", wintypes.DWORD),
                        ("ActiveProcesses", wintypes.DWORD), ("TotalTerminatedProcesses", wintypes.DWORD)]

        class ThreadEntry(ctypes.Structure):            # THREADENTRY32
            _fields_ = [("dwSize", wintypes.DWORD), ("cntUsage", wintypes.DWORD), ("th32ThreadID", wintypes.DWORD),
                        ("th32OwnerProcessID", wintypes.DWORD), ("tpBasePri", wintypes.LONG), ("tpDeltaPri", wintypes.LONG),
                        ("dwFlags", wintypes.DWORD)]

        k.CreateJobObjectW.argtypes, k.CreateJobObjectW.restype = [ctypes.c_void_p, wintypes.LPCWSTR], wintypes.HANDLE
        k.SetInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD]
        k.SetInformationJobObject.restype = wintypes.BOOL
        k.QueryInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD, ctypes.c_void_p]
        k.QueryInformationJobObject.restype = wintypes.BOOL
        k.AssignProcessToJobObject.argtypes, k.AssignProcessToJobObject.restype = [wintypes.HANDLE, wintypes.HANDLE], wintypes.BOOL
        k.TerminateJobObject.argtypes, k.TerminateJobObject.restype = [wintypes.HANDLE, wintypes.UINT], wintypes.BOOL
        k.OpenProcess.argtypes, k.OpenProcess.restype = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD], wintypes.HANDLE
        k.CloseHandle.argtypes, k.CloseHandle.restype = [wintypes.HANDLE], wintypes.BOOL
        k.CreateToolhelp32Snapshot.argtypes = [wintypes.DWORD, wintypes.DWORD]
        k.CreateToolhelp32Snapshot.restype = wintypes.HANDLE
        k.Thread32First.argtypes = k.Thread32Next.argtypes = [wintypes.HANDLE, ctypes.c_void_p]
        k.Thread32First.restype = k.Thread32Next.restype = wintypes.BOOL
        k.OpenThread.argtypes, k.OpenThread.restype = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD], wintypes.HANDLE
        k.ResumeThread.argtypes, k.ResumeThread.restype = [wintypes.HANDLE], wintypes.DWORD
        _WIN = types.SimpleNamespace(ctypes=ctypes, k=k, ExtendedLimit=ExtendedLimit, Accounting=Accounting, ThreadEntry=ThreadEntry)
    return _WIN


def _win_error(what: str) -> RuntimeError:
    """RuntimeError that says what could not be done and why (Windows' own words for the error of the call that just failed)."""
    w = _win()
    return RuntimeError(f"cannot {what}: {w.ctypes.WinError(w.ctypes.get_last_error())}")


class _Job:
    """A Windows job object that ends every process in it when the last handle to it is closed, and holds every process in it at one
    priority class (module docstring, "run_solver").

    The sequence run_solver uses: job = _Job(class); start the program with _CREATE_SUSPENDED; job.assign(pid); job.resume(pid); ...
    job.active() to see who is left; job.end() to end them; job.close() in a finally. self.handle is the only handle to the job and
    is not inheritable, so no child holds it: it closes when this process closes it or dies, and then the system ends the members.

    The class: JOB_OBJECT_LIMIT_PRIORITY_CLASS (0x20) with PriorityClass = `priority_class`, the class run_solver was asked for
    (BELOW_NORMAL_PRIORITY_CLASS 0x4000, or IDLE_PRIORITY_CLASS 0x40: _windows_class). "Causes all processes associated with the job
    to use the same priority class": a member that calls SetPriorityClass (the call succeeds and changes nothing), and a member that
    starts a child with another class in its creation flags (the child is born at the job's), cannot leave it. Without the limit a
    program that was started below normal could raise itself or its children; MATLAB's graphics hosts did. Both hold for a program that
    is not elevated: tests/test_femm_queue.py shows it with a stand-in program, and its control without the job.

    `who` is what the error messages call the process that goes into the job: "the program" for run_solver's program (the default),
    "this worker" for a pool worker that joins a job of its own (_join_own_job), whose boot-failed text carries them.
    """

    def __init__(self, priority_class: int, who: str = "the program"):
        w = _win()
        self.who = who
        self.handle = w.k.CreateJobObjectW(None, None)      # no security attributes: the handle is not inheritable
        if not self.handle:
            raise _win_error(f"create a job object for {who}")
        limits = w.ExtendedLimit()
        limits.BasicLimitInformation.LimitFlags = 0x2000 | 0x20     # JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE | JOB_OBJECT_LIMIT_PRIORITY_CLASS
        limits.BasicLimitInformation.PriorityClass = priority_class
        if not w.k.SetInformationJobObject(self.handle, 9, w.ctypes.byref(limits), w.ctypes.sizeof(limits)):   # 9 = extended limits
            problem = _win_error("make the job end its processes when it is closed and hold them at their priority class")
            self.close()
            raise problem

    def assign(self, pid: int) -> None:
        """Put the process `pid` into the job: run_solver's program (our own suspended child: its Popen holds a handle, so the number is still
        its) or, for a pool worker, this very process."""
        w = _win()
        process = w.k.OpenProcess(0x0100 | 0x0001, False, pid)      # PROCESS_SET_QUOTA | PROCESS_TERMINATE: what assigning needs
        if not process:
            raise _win_error(f"open {self.who} to put it into its job")
        try:
            if not w.k.AssignProcessToJobObject(self.handle, process):
                problem = _win_error(f"put {self.who} into its job")
                raise RuntimeError(f"{problem} (a job that already holds it may allow no other)")
        finally:
            w.k.CloseHandle(process)

    def resume(self, pid: int) -> None:
        """Let the suspended program run: find its one thread with a Toolhelp snapshot and ResumeThread it.

        Popen closes the thread handle that CreateProcess returns, so the thread has to be found again. The documented way is
        a snapshot of the system's threads, filtered by owner. (NtResumeProcess would do it in one call, but it is not documented.)
        """
        w = _win()
        snapshot = w.k.CreateToolhelp32Snapshot(0x4, 0)             # TH32CS_SNAPTHREAD
        if snapshot is None or snapshot == w.ctypes.c_void_p(-1).value:
            raise _win_error("list the threads to start the program")
        counts: list[int] = []                                      # the suspend count each ResumeThread returned
        try:
            entry = w.ThreadEntry()
            entry.dwSize = w.ctypes.sizeof(entry)
            more = w.k.Thread32First(snapshot, w.ctypes.byref(entry))
            while more:
                if entry.th32OwnerProcessID == pid:
                    thread = w.k.OpenThread(0x0002, False, entry.th32ThreadID)      # THREAD_SUSPEND_RESUME
                    if not thread:
                        raise _win_error("open the program's thread to start it")
                    try:
                        counts.append(w.k.ResumeThread(thread))
                    finally:
                        w.k.CloseHandle(thread)
                more = w.k.Thread32Next(snapshot, w.ctypes.byref(entry))
        finally:
            w.k.CloseHandle(snapshot)
        if counts != [1]:                                           # a new suspended process has exactly one thread, suspended once
            raise RuntimeError(f"cannot start the program: expected to resume one thread that was suspended once, resumed {counts}")

    def active(self) -> int:
        """How many processes of the job are still running."""
        w = _win()
        info = w.Accounting()
        if not w.k.QueryInformationJobObject(self.handle, 1, w.ctypes.byref(info), w.ctypes.sizeof(info), None):   # 1 = accounting
            raise _win_error("read how many processes of the job are still running")
        return info.ActiveProcesses

    def end(self, wait_s: float = 10.0) -> None:
        """End every process in the job, orphans included (TerminateJobObject), and wait until none is left."""
        w = _win()
        if not w.k.TerminateJobObject(self.handle, 1):
            raise _win_error("end the program's job")
        deadline = time.monotonic() + wait_s
        while self.active() > 0:
            if time.monotonic() >= deadline:
                raise RuntimeError(f"{self.active()} process(es) of the program are still running {wait_s:g} s after its job was ended")
            time.sleep(0.02)

    def close(self) -> None:
        """Close the handle. With processes still in the job this ends them (that is what the limit is for). Safe to call twice."""
        if self.handle:
            _win().k.CloseHandle(self.handle)
            self.handle = None


def _end_unstarted(proc: subprocess.Popen) -> None:
    """The program was created suspended and could not be put into its job or started: end it. Nothing has run, nothing is below it."""
    proc.kill()
    _reap(proc)


def _run_in_job(args: list, cwd, env: dict, priority: str, timeout_s: float) -> SolverRun:
    """run_solver on Windows: the program, and everything it starts, lives in a job that dies with this process and holds every process in it at the class asked for."""
    job = _Job(_windows_class(priority))                        # before the program exists: if this fails, nothing has been started
    try:
        started = time.monotonic()
        proc = None                                             # set inside the try: a Ctrl-C right after Popen cannot lose the program
        try:
            proc = subprocess.Popen(args, cwd=cwd, env=env, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    encoding="utf-8", errors="replace", creationflags=_windows_class(priority) | _CREATE_SUSPENDED)
            job.assign(proc.pid)                                # while it is suspended: nothing it starts can be born outside
            job.resume(proc.pid)
        except BaseException:
            if proc is not None:                                # (None: Popen itself failed, so there is no program to end)
                _end_unstarted(proc)
            raise
        deadline = started + timeout_s
        try:
            out, err = _communicate(proc, deadline, timeout_s)  # until the program's own process has ended ...
            while job.active() > 0:                             # ... and then until what it started has ended too
                if time.monotonic() >= deadline:
                    raise subprocess.TimeoutExpired(args, timeout_s)
                time.sleep(0.05)
            return SolverRun(proc.returncode, out, err, time.monotonic() - started)
        except subprocess.TimeoutExpired:
            job.end()                                           # the program and every orphan below it
            out, err = _last_output(proc)
            return SolverRun(None, out, err, time.monotonic() - started)
        except BaseException:                                   # Ctrl-C, or an error here: never leave the program running
            job.end()
            _reap(proc)
            raise
    finally:
        job.close()                                             # anything still in the job is ended by this, whatever happened above


def _run_in_group(args: list, cwd, env: dict, priority: str, timeout_s: float) -> SolverRun:
    """run_solver off Windows: the program leads its own process group, which a timeout or Ctrl-C ends. A caller that is itself
    ended without clean-up can leave the program running: there is no job to close (module docstring)."""
    started = time.monotonic()
    proc = subprocess.Popen(args, cwd=cwd, env=env, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            encoding="utf-8", errors="replace", **_popen_options(priority))
    try:
        out, err = _communicate(proc, started + timeout_s, timeout_s)
        return SolverRun(proc.returncode, out, err, time.monotonic() - started)
    except subprocess.TimeoutExpired:
        _end_program(proc)
        out, err = _last_output(proc)
        return SolverRun(None, out, err, time.monotonic() - started)
    except BaseException:                                       # Ctrl-C, or an error here: never leave the program running
        _end_program(proc)
        raise


def run_solver(cmd, *, label: str, cwd, timeout_s: float, threads: int | None = None,
               wait: bool = True, poll_s: float = 30.0, lock_dir=None) -> SolverRun:
    """Run one solver program under the solver slot, at below-normal priority, for at most timeout_s seconds. On Windows every process of
    the program's job is held at that class (_Job), whatever it asks for itself.

    cmd       a list: [program, argument, ...] (never one string: a path with a space in it would be cut in two)
    label     what a waiting run prints as the holder, e.g. "bench B4 ccx"
    cwd       the folder the program runs in (a job folder of its own, so two programs never share file names)
    timeout_s seconds the PROGRAM and everything it starts may run, counted from its start (the wait for the slot does not count).
              Required: every run is under a timeout (R10). When it is over, all of it is ended: through the program's job on
              Windows (orphans included), through its process group elsewhere.
    threads   how many threads the program may use, 1 to [femm] workers (8). It sets OMP_NUM_THREADS and CalculiX's
              CCX_NPROC_EQUATION_SOLVER and CCX_NPROC_STIFFNESS to that number, in the program's environment only. None means
              the maximum, so a caller that does not choose is still capped.
    wait, poll_s, lock_dir   as for solver_slot()
    Returns SolverRun(returncode, stdout, stderr, seconds). returncode is the program's own exit code, and `seconds` runs until the
    whole tree has ended (on Windows; see the module docstring). If it was ended on timeout, returncode is None (.timed_out is True)
    and stdout and stderr hold what it had written until then. Raises ValueError or TypeError, before anything is taken or started,
    for a bad argument; RuntimeError if the slot is already held by this thread or by the run that started this process, or (Windows)
    if the program cannot be put into its job: the program is then ended before it has run anything.
    The limit off Windows: the program leads its own process group, which a timeout or a Ctrl-C ends, but a caller that is itself
    ended without clean-up (kill -9) can leave the program running under a free slot. On Windows its job makes that impossible.
    """
    if isinstance(cmd, (str, bytes, os.PathLike)) or len(cmd) == 0:
        raise TypeError(f"cmd must be a list [program, argument, ...], not {cmd!r}: one string would be split at every space")
    args = [os.fspath(part) for part in cmd]
    count = _limit("threads", threads)
    if isinstance(timeout_s, bool) or not isinstance(timeout_s, (int, float)) or not timeout_s > 0:
        raise ValueError(f"timeout_s={timeout_s!r}: a positive number of seconds (rule R10: every run is under a timeout)")
    priority = _priority_name()
    with solver_slot(label, wait=wait, poll_s=poll_s, lock_dir=lock_dir):
        env = dict(os.environ)                                  # a copy: this process's own environment is not touched
        env.update(OMP_NUM_THREADS=str(count), CCX_NPROC_EQUATION_SOLVER=str(count), CCX_NPROC_STIFFNESS=str(count))
        if sys.platform == "win32":
            return _run_in_job(args, cwd, env, priority, timeout_s)
        return _run_in_group(args, cwd, env, priority, timeout_s)


# ============================================================================================== the pool: a worker's side
class WorkerTraceback(Exception):
    """The traceback of an error raised in a worker, as text. It is the __cause__ of the exception that map() raises."""

    def __init__(self, text: str):
        super().__init__(text)
        self.text = text

    def __str__(self) -> str:
        return "\n" + self.text


def _describe(error: BaseException) -> tuple:
    """An exception in a form that crosses a pipe: (pickled exception or None, its type's name, its message, its traceback)."""
    try:
        blob = pickle.dumps(error)
    except Exception:                                           # some exceptions cannot be pickled; the text still can
        blob = None
    return blob, type(error).__name__, str(error), "".join(traceback.format_exception(error))


def _rebuild(description: tuple) -> BaseException:
    """The exception that _describe() packed, with the worker's traceback as its __cause__."""
    blob, name, text, trace = description
    error = None
    if blob is not None:
        try:
            error = pickle.loads(blob)
        except Exception:
            error = None
    if error is None:
        error = RuntimeError(f"{name}: {text}")
    error.__cause__ = WorkerTraceback(trace)
    return error


_OWN_JOB = None                                 # Windows, in a pool worker: the worker's job. Its only handle, kept for the worker's whole life


def _join_own_job(priority: str) -> None:
    """Windows, in a pool worker, before it starts any process: put THIS process into a job of its own, of the kind run_solver's programs live in
    (_Job: kill every member when the last handle is closed, hold every member at the pool's priority class), and keep the job's only handle in
    _OWN_JOB for the worker's whole life. Closing it would end the worker itself. When the worker ends, however it ends (the pool ending it with
    Process.kill, TerminateProcess, End task, a crash, its own watch thread), the system closes the handle and ends every process the worker started, FEMM and its
    fkn.exe included: a member's children are born members, and nothing here asks to break away from the job (femm_session starts FEMM with a plain Popen).
    The handle is made with no security attributes (_Job), so it is not inheritable and no child holds the job open.

    Why a job and not a walk of the worker's tree (taskkill /T, which this replaced on Windows, in the watch thread and in the pool's timeout kill): a
    walk finds what exists when it looks. A FEMM that the worker's main thread created a moment later, before the worker was gone, was not in it and
    outlived its worker (the timeout trial of B9, 2026-11-06). A job has no moment to be early for, and it holds what the worker started and nothing
    else, where a walk goes by the parent numbers the processes recorded.

    RuntimeError, naming the cause, if the job cannot be made or this process cannot be put into it (the caller's own job may allow no other): the
    worker then tells its owner ("boot-failed") and starts nothing, not FEMM and not the caller's initializer. A worker is never run unprotected, silently:
    the rule of run_solver."""
    global _OWN_JOB
    if _OWN_JOB is not None:
        raise RuntimeError("this process is already in a job of its own")
    try:
        job = _Job(_windows_class(priority), who="this worker")
        try:
            job.assign(os.getpid())
        except Exception:
            job.close()                                         # nothing is in it, so closing it ends nobody
            raise
    except Exception as error:
        raise RuntimeError(f"a pool worker cannot join a job of its own, so it does not start FEMM or the caller's initializer (a worker is never run "
                           f"unprotected, silently): {error}") from error
    _OWN_JOB = job


def _watch_owner() -> None:
    """Start a thread that ends this worker the moment the process that started it is gone.

    Without it a killed owner would leave its workers (each with a FEMM) computing after the slot is free, and the next run
    would start on top of them. multiprocessing.parent_process() is a handle on the owner that can be waited on.

    Windows: the thread ends the worker at once (os._exit) and walks no tree. Whatever the worker started is in the worker's own job
    (_join_own_job), and the system ends it with the worker. Elsewhere the worker leads its own process group (setsid), and the thread
    ends the group (_end_tree: killpg), which takes everything below the worker with it.
    """
    owner = multiprocessing.parent_process()
    if owner is None:
        return

    def watch() -> None:
        owner.join()                                            # returns when the owner process ends, however it ends
        try:
            if sys.platform != "win32":
                _end_tree(os.getpid())                          # this worker and what is below it (it ends this thread too)
        finally:
            os._exit(1)                                         # (on Windows this is all: the job ends the rest)

    threading.Thread(target=watch, name="watch-owner", daemon=True).start()


def _worker_main(conn, priority: str, use_femm: bool, initializer, initargs) -> None:
    """What every pool worker runs. Messages to its owner: ("ready", pid), ("boot-failed", traceback), ("done", index, result),
    ("error", index, description). From the owner: (index, pickled fn, pickled item), or None for "finish up and quit". The function
    and the item travel as bytes and are unpickled in here, inside the try: a function that this worker cannot import (a module
    that exists only in the owner) is then an error with its cause, not a worker that dies for no visible reason.

    On Windows the worker joins a job of its own (_join_own_job) before it starts anything: if it cannot, it says "boot-failed" with the
    cause and calls neither femm_session.open_femm() nor the caller's initializer."""
    if sys.platform != "win32":
        os.setsid()                                             # lead our own process group, so _end_tree() can end it with FEMM
    signal.signal(signal.SIGINT, signal.SIG_IGN)                # a Ctrl-C in the terminal goes to the owner, which ends us
    _lower_self(priority)                                       # first thing: nothing this worker starts is ever above it
    _watch_owner()
    femm_open = False
    try:
        if sys.platform == "win32":
            _join_own_job(priority)                             # before anything is started: FEMM, fkn.exe and whatever the caller's code starts are born in the job
        if use_femm:
            import femm_session
            femm_session.open_femm()                            # lowers the priority first, then starts this worker's own FEMM
            femm_open = True
        if initializer is not None:
            initializer(*initargs)                              # the caller's: after FEMM is open
        conn.send(("ready", os.getpid()))
    except BaseException:
        try:
            conn.send(("boot-failed", traceback.format_exc()))
        except (OSError, ValueError):
            pass
        _close_femm(femm_open)
        return
    try:
        while True:
            try:
                job = conn.recv()
            except (EOFError, OSError):                         # the owner is gone
                break
            if job is None:                                     # "finish up and quit"
                break
            index, fn_bytes, item_bytes = job
            try:
                fn, item = pickle.loads(fn_bytes), pickle.loads(item_bytes)
                conn.send(("done", index, fn(item)))
            except BaseException as error:                      # SystemExit too: whatever it was, the owner is told
                conn.send(("error", index, _describe(error)))
    finally:
        _close_femm(femm_open)


def _close_femm(opened: bool) -> None:
    if opened:
        import femm_session
        try:
            femm_session.close_femm()
        except Exception:
            pass


# ============================================================================================== the pool: the owner's side
class _Worker:
    """One worker as its owner sees it: the process, the pipe to it, and what it is doing."""

    def __init__(self, process, conn):
        self.process, self.conn = process, conn
        self.state = "starting"                                 # starting -> idle <-> busy; dead is removed from the pool at once
        self.task = None                                        # (index, failures so far) while busy
        self.deadline = None                                    # when the task in hand is over time (time.monotonic), or None
        self.born = time.monotonic()


class Pool:
    """The workers of one femm_pool: what `with femm_queue.femm_pool(...) as pool` provides. Use .map(); .pids and .workers say who.

    Not thread-safe: one thread calls map(), one call at a time. femm_pool() makes it, starts it and closes it.
    """

    def __init__(self, label, workers, task_timeout_s, start_timeout_s, initializer, initargs, use_femm, priority):
        self.label, self.workers = label, workers
        self.task_timeout_s, self._start_timeout_s = task_timeout_s, start_timeout_s
        self._initializer, self._initargs, self._use_femm, self._priority = initializer, initargs, use_femm, priority
        self._workers: list[_Worker] = []
        self._closing = False
        self._pending: deque = deque()
        self._results: list = []
        self._remaining = 0

    @property
    def pids(self) -> list[int]:
        """The process ids of the workers that are running now."""
        return [w.process.pid for w in self._workers]

    def _say(self, message: str) -> None:
        print(f"[femm_pool {self.label}] {message}", flush=True)

    # ---------------------------------------------------------------- starting and ending workers
    def _spawn(self) -> _Worker:
        """Start one worker. It is "starting" until it says ready."""
        owner_end, worker_end = _CTX.Pipe(duplex=True)
        process = _CTX.Process(target=_worker_main, name=f"femm_pool {self.label}",
                               args=(worker_end, self._priority, self._use_femm, self._initializer, self._initargs))
        process.start()
        worker_end.close()                  # keep only our end: then a worker that dies shows up here as end-of-file
        worker = _Worker(process, owner_end)
        self._workers.append(worker)
        return worker

    def _end_worker(self, worker: _Worker) -> None:
        """End this worker and everything it started (its FEMM and the fkn.exe of each solve, a solver of its own).

        Windows: Process.kill(), which is TerminateProcess through the handle the pool holds on the worker, never through its number, and the
        worker's own job (_join_own_job) ends the rest. No taskkill: its /T walk goes by the parent numbers the processes recorded and could
        take an orphan whose dead parent's number the worker now holds (the module docstring has the whole argument). Elsewhere: SIGKILL to the
        process group the worker leads (_end_tree)."""
        process = worker.process
        assert process.pid != os.getpid(), "a pool never ends its own owner"
        said = ""
        for _ in range(2):
            if not process.is_alive():
                return
            if sys.platform == "win32":
                process.kill()
            else:
                code, text = _end_tree(process.pid)
                said = f"; the signal said {code}: {text}"
            process.join(10)
        if process.is_alive():
            raise RuntimeError(f"femm_pool {self.label!r}: could not end worker pid {process.pid}{said}")

    def _discard(self, worker: _Worker) -> None:
        """Forget a worker that is not running: take it out of the pool and release its handles."""
        worker.state = "dead"
        if worker in self._workers:
            self._workers.remove(worker)
        try:
            worker.conn.close()
            worker.process.close()
        except (OSError, ValueError):
            pass

    def _retire(self, worker: _Worker) -> None:
        self._end_worker(worker)
        self._discard(worker)

    def _start(self) -> None:
        """Start `workers` workers and wait until every one is ready. RuntimeError if one fails or does not start in time."""
        for _ in range(self.workers):
            self._spawn()
        self._say(f"starting {self.workers} worker(s), pids {self.pids}, at {self._priority} priority"
                  + (", each with its own FEMM" if self._use_femm else ""))
        while any(w.state == "starting" for w in self._workers):
            self._check_start_deadlines()
            self._pump(0.5)
        self._say("all workers are ready")

    def _check_start_deadlines(self) -> None:
        now = time.monotonic()
        late = [w for w in self._workers if w.state == "starting" and now - w.born > self._start_timeout_s]
        if late:
            raise RuntimeError(
                f"femm_pool {self.label!r}: worker(s) {[w.process.pid for w in late]} did not start within "
                f"{self._start_timeout_s:g} s. A FEMM that never answers is the usual cause. They are ended.")

    # ---------------------------------------------------------------- news from the workers
    def _pump(self, timeout: float) -> None:
        """Wait up to `timeout` seconds for news from the workers (a result, a message, a death) and act on it."""
        by_handle = {}                                          # what to wait on -> its worker: the pipe, and the process itself
        for w in self._workers:                                 # (a process that dies is seen even if something below it kept the pipe open)
            by_handle[w.conn] = by_handle[w.process.sentinel] = w
        if not by_handle:                                       # cannot happen while a map() runs; never spin if it does
            time.sleep(timeout)
            return
        touched: list[_Worker] = []
        for handle in multiprocessing.connection.wait(list(by_handle), timeout):
            if by_handle[handle] not in touched:
                touched.append(by_handle[handle])
        for worker in touched:
            if worker.state == "dead":
                continue
            heard, closed = [], False
            try:
                while worker.conn.poll(0):                      # everything it has said ...
                    heard.append(worker.conn.recv())
            except (EOFError, OSError):                         # ... and then the pipe is closed: the worker ended
                closed = True
            # Only the reading is inside the try. _handle() raises the task's own exception, which may well be an OSError (a model
            # file that is not there): that is the task's error, and must not be mistaken for a worker that died.
            for message in heard:
                self._handle(worker, message)
            if closed or not worker.process.is_alive():
                self._lost(worker)

    def _handle(self, worker: _Worker, message: tuple) -> None:
        kind = message[0]
        if kind == "ready":
            worker.state = "idle"
        elif kind == "boot-failed":
            raise RuntimeError(f"femm_pool {self.label!r}: a worker failed while starting (opening its FEMM, or the caller's initializer):\n"
                               f"{message[1]}")
        elif kind == "done":
            self._results[message[1]] = message[2]
            self._remaining -= 1
            worker.state, worker.task, worker.deadline = "idle", None, None
        elif kind == "error":
            worker.state, worker.task, worker.deadline = "idle", None, None     # the worker is fine: its task raised, that is all
            raise _rebuild(message[2])

    def _lost(self, worker: _Worker) -> None:
        """A worker ended by itself (not because we ended it). Its task, if it had one, counts as a failure."""
        was, task = worker.state, worker.task
        worker.process.join(5)
        if worker.process.is_alive():                           # it closed its pipe but is still there: it is ours, end it
            self._end_worker(worker)
        code = worker.process.exitcode
        pid = worker.process.pid
        self._discard(worker)
        if was == "starting":
            raise RuntimeError(f"femm_pool {self.label!r}: worker pid {pid} died while starting (exit code {code}); its own traceback, "
                               "if it had one, is printed above. A script without the `if __name__ == \"__main__\":` guard does this: "
                               "every worker re-runs its top level.")
        self._say(f"worker pid {pid} ended by itself (exit code {code}); starting a fresh one")
        if not self._closing:
            self._spawn()
        if was == "busy":
            self._failed(*task, {"status": "died"}, f"its worker died (exit code {code})")

    def _failed(self, index: int, failures: int, final: dict, why: str) -> None:
        """An attempt at an item failed. The first failure gets one more try; the second ends the item with `final`."""
        if failures == 0:
            self._pending.appendleft((index, 1))
            self._say(f"item {index}: {why}; trying it once more")
        else:
            self._results[index] = final
            self._remaining -= 1
            self._say(f"item {index}: {why} again; giving up on it, its result is {final}")

    def _timed_out(self, worker: _Worker) -> None:
        index, failures = worker.task
        self._say(f"item {index} ran over {self.task_timeout_s:g} s in worker pid {worker.process.pid}: "
                  "ending that worker and everything below it")
        self._retire(worker)
        self._spawn()                                           # a fresh worker takes its place (the task may go to any free one)
        self._failed(index, failures, {"status": "timeout"}, f"over {self.task_timeout_s:g} s")

    # ---------------------------------------------------------------- map
    def map(self, fn, items) -> list:
        """fn(item) for every item, on the workers, in parallel. Returns the results in the order of the items.

        fn must be a module-level function and the items picklable (a worker is a new Python that imports fn by name; a lambda
        raises TypeError here, and a function a worker cannot import is raised here as the worker's ImportError).
        An item whose task is over task_timeout_s twice has the result {"status": "timeout"}; one whose worker died twice,
        {"status": "died"}. An exception in a task is raised here, with the worker's traceback as its __cause__, and the
        tasks still running are ended. Can be called more than once.
        """
        items = list(items)
        if self._closing:
            raise RuntimeError("this pool is closed")
        if not items:
            return []
        try:
            fn_bytes = pickle.dumps(fn)                         # by name: a worker imports the module and finds the function
        except Exception as error:
            raise TypeError(f"fn must be a function defined at the top level of a module that a new Python can import "
                            f"(a lambda or a function inside a function cannot be sent to a worker): {error}") from error
        while len(self._workers) < self.workers:                # a worker was lost in an earlier, aborted call
            self._spawn()
        self._results = [None] * len(items)
        self._pending = deque((index, 0) for index in range(len(items)))
        self._remaining = len(items)
        try:
            while self._remaining:
                self._dispatch(fn_bytes, items)
                self._pump(self._wait_time())
                self._check_deadlines()
        except BaseException:
            self._abandon()
            raise
        return self._results

    def _dispatch(self, fn_bytes: bytes, items: list) -> None:
        """Give every idle worker a task, while there are tasks."""
        for worker in list(self._workers):
            if worker.state != "idle" or not self._pending:
                continue
            index, failures = self._pending.popleft()
            item_bytes = pickle.dumps(items[index])             # an item that cannot be pickled raises here, in the caller
            try:
                worker.conn.send((index, fn_bytes, item_bytes))
            except (OSError, EOFError):                             # the worker died while idle: not this task's fault
                self._pending.appendleft((index, failures))
                self._lost(worker)
                continue
            worker.state, worker.task = "busy", (index, failures)
            worker.deadline = None if self.task_timeout_s is None else time.monotonic() + self.task_timeout_s

    def _wait_time(self) -> float:
        """How long to wait for news: until the nearest deadline, at most a second (so a Ctrl-C is never kept waiting)."""
        deadlines = [w.deadline for w in self._workers if w.state == "busy" and w.deadline is not None]
        if not deadlines:
            return 1.0
        return max(0.0, min(min(deadlines) - time.monotonic(), 1.0))

    def _check_deadlines(self) -> None:
        self._check_start_deadlines()
        now = time.monotonic()
        for worker in list(self._workers):
            if worker.state == "busy" and worker.deadline is not None and now >= worker.deadline:
                self._timed_out(worker)

    def _abandon(self) -> None:
        """A map() is being left by an exception: end the workers that are still busy (their results are not wanted)."""
        for worker in list(self._workers):
            if worker.state == "busy":
                self._retire(worker)

    # ---------------------------------------------------------------- closing
    def close(self) -> None:
        """Finish: idle workers are told to quit (they close their FEMM), any other worker is ended, nothing is left running."""
        self._closing = True
        polite = [w for w in self._workers if w.state == "idle"]
        rough = [w for w in self._workers if w.state != "idle"]
        for worker in polite:
            try:
                worker.conn.send(None)
            except (OSError, ValueError):
                pass
        problems: list[Exception] = []
        for worker in rough:
            try:
                self._end_worker(worker)
            except RuntimeError as error:
                problems.append(error)
        deadline = time.monotonic() + 30
        for worker in polite:
            worker.process.join(max(0.0, deadline - time.monotonic()))
            if worker.process.is_alive():
                self._say(f"worker pid {worker.process.pid} did not quit in time: ending it")
                try:
                    self._end_worker(worker)
                except RuntimeError as error:
                    problems.append(error)
        for worker in list(self._workers):
            self._discard(worker)
        if problems:
            raise problems[0]


@contextmanager
def femm_pool(label: str, workers: int | None = None, task_timeout_s: float | None = None, initializer=None, initargs=(),
              *, femm: bool = True, start_timeout_s: float = _START_TIMEOUT_S, wait: bool = True, poll_s: float = 30.0,
              lock_dir=None):
    """Take the solver slot, start the workers, yield a Pool (use pool.map(fn, items) -> list), and on the way out end every worker.

    workers          how many worker processes: 1 to config.toml's [femm] workers (8); None means that maximum. More is refused.
    task_timeout_s   a task over this many seconds (counted from the moment a ready worker is handed it) has its worker's process
                     tree ended and is tried once more; the second time its result is {"status": "timeout"}. None means no limit
                     (R10 wants one: pass it).
    initializer      called as initializer(*initargs) in every worker once, AFTER its FEMM is open (a module-level function)
    femm             True: each worker first opens its own FEMM with femm_session.open_femm() and closes it at its end.
                     False: no FEMM (the tests, and sweeps whose tasks run something else); the workers still run below normal.
    start_timeout_s  a worker that is not ready after this long did not start: the pool is not opened, RuntimeError. This is the net
                     for a FEMM that never answers (femm_session.open_femm() has no timeout of its own)
    wait, poll_s, lock_dir   as for solver_slot()
    All refusals (a bad number, a priority that is not below normal) come before the slot is taken or anything is started.
    """
    count = _limit("workers", workers)
    priority = _priority_name()
    if task_timeout_s is not None and (isinstance(task_timeout_s, bool) or not isinstance(task_timeout_s, (int, float))
                                       or not task_timeout_s > 0):
        raise ValueError(f"task_timeout_s={task_timeout_s!r}: None, or a positive number of seconds")
    with solver_slot(label, wait=wait, poll_s=poll_s, lock_dir=lock_dir):
        pool = Pool(label, count, task_timeout_s, start_timeout_s, initializer, tuple(initargs), femm, priority)
        try:
            pool._start()
            yield pool
        finally:
            pool.close()


# ============================================================================================== the command line
def main(argv: list[str] | None = None) -> int:
    """`python femm_queue.py [lock folder]`: say who holds the solver slot. Exit status 0 if it is free, 1 if it is held."""
    args = sys.argv[1:] if argv is None else argv
    info = holder(args[0] if args else None)
    if info is None:
        print("the solver slot is free")
        return 0
    print(f"the solver slot is held by {_who(info)}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
