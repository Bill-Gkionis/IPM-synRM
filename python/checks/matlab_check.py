"""matlab_check: `matlab -batch` loads a copy of the motor's Simulink plant, simulates a small model it builds, and we see what it leaves.

What it does, and why
    MATLAB and Simulink are where the control and the inverter are modelled (the plant, Improved_plant.slx, is the motor as the
    controller sees it). This row proves that a script can start MATLAB, use Simulink, get an answer back and let go, without a
    person at the keyboard: `matlab -batch "statement"` runs one statement without a window and exits with a code.
    The script it runs (written into the case folder by this module, check_case.m) does two things:
      1. the plant: it loads motor structure.mat into the workspace (the parameters the model needs), loads a COPY of
         Improved_plant.slx with load_system, does NOT simulate or update it, and closes it without saving;
      2. a small model: it builds Sine -> Gain(3) -> To Workspace, simulates one second at 0.01 s, and checks that the output is
         exactly 3 times the input and that the input is 2 sin(5 t) + 0.5 (to 1e-12).
    Any error in the script makes MATLAB exit with a nonzero code. The row passes when the exit code is 0, the script reached its
    last line (it writes OK=1 into result.txt), the plant copy and the two reference files are byte for byte what they were, the
    job (below) emptied, and MATLAB computed with at most 8 threads (below).
    The threads (R10: at most 8 workers or solver threads): run_solver's threads=4 sets OMP_NUM_THREADS, which MathWorks does not
    document as limiting MATLAB's computational threads; their default is the number of physical cores (12 on this PC). So MATLAB is started with
    -singleCompThread (MathWorks: "Limit MATLAB to a single computational thread"), as the bench's B8 does, and the script writes
    what maxNumCompThreads says into result.txt (MAXNUMCOMPTHREADS=1). The row fails if that is missing, below 1 or above [femm] workers.

The two reference files (config.toml [paths] simulink_plant, motor_structure_mat) live in P\\reference, committed. They are
hashed before and after, and MATLAB only ever sees copies in the case folder, so no run can change the originals.

What the job held (R10, femm_queue.run_solver)
    run_solver puts the program into a Windows job and returns only when the job is EMPTY, not when the program's own process
    exits: matlab.exe is a launcher that starts MATLAB.exe, and the slot must not be free while MATLAB computes. The price is that a
    helper process that stays running after MATLAB would keep run_solver waiting until the timeout. The design decision is for it to
    be MEASURED here, not worked around: while MATLAB runs, a thread takes a list of the PC's processes twenty times a second (a read,
    like Task Manager's), follows every process that descends from the program run_solver started, and records when each was first
    and last seen. A process is known by its pid AND its creation time (read through a handle opened with
    PROCESS_QUERY_LIMITED_INFORMATION, which can read a process and cannot end or change it; this file asks for no other right),
    because Windows hands the pid of an ended process to a new one. If the list cannot be taken, the summary says how many looks
    failed, and never reads as a job that held nothing. The summary says which processes joined the job and how long after the
    last matlab.exe the job was empty.
    Only the "warm" case can be measured here: MathWorksServiceHost.exe starts at logon (it is a Run-key item of this user), so it is
    already running OUTSIDE any job when MATLAB starts; nothing here ends it or any process this run did not start.
    The summary goes to the measured text and, whole, to job_members.txt in the case folder.

Numbers to expect (2026-11-06, MATLAB R2025b Update 5, warm; four runs, the same each time)
    * the row takes 11 to 12 s; exit code 0; the plant has 2 parameter variables in motor structure.mat and a stop time of 0.5 s
      (untouched); the gain model gives 101 samples, a gain error of 0 and a sine error of 7.6e-15; MATLAB reports 1 computational
      thread (the row has been run with -singleCompThread since the final review of the tools plan);
    * load_system prints a warning that the target file veristand.tlc cannot be found (the plant was saved with NI VeriStand's
      code-generation target, which is not installed; it matters for building code, not for loading the model);
    * the job held 3 to 5 processes, all of MATLAB's own: matlab.exe (the launcher: it stays until MATLAB.exe is done, so its
      exit code is MATLAB's), MATLAB.exe, and MathWorksServiceHost.exe processes that MATLAB starts as ITS children (one that lasts
      until about 0.5 s before MATLAB.exe ends, two that last 0.1 to 0.2 s). The MathWorksServiceHost that runs all day (started at
      logon) is not among them. The job was empty less than 0.05 s (one look) after the last matlab.exe: nothing stayed behind.
    * a process that lives less than the 0.05 s between two looks can be missed; the list is "at least these".

Self-check (the answers follow).
  1. Why is the exit code not enough to say the script ran to its end?
  2. Why does the plant part not simulate the plant?
  3. What would be seen here if MathWorksServiceHost.exe had been started by MATLAB instead of at logon?
Answers.
  1. The exit code is the launcher's. Here matlab.exe waits for MATLAB.exe and passes its status on (the job tracker saw it stay
     for the whole run), but that is a property of this launcher; a launcher that returned at once would report 0 for a script that
     had not started. The script's last line writes OK=1 into a file: only a run that got there has it.
  2. Simulating the plant takes minutes and needs every parameter of the motor set; the brief asks only for the model to load
     (which proves the file, the library links and the release agree). Task 9's check simulated it for 2 ms.
  3. That is the cold case, which this PC cannot show (the host starts at logon). A long-lived host that MATLAB had to start would
     become a member of the job and would not end with MATLAB (that is what makes it a service), so run_solver would wait until the
     timeout (5 minutes), end it with the job and report timed out; this row would fail and name the process. It would be a case to
     report, not to hide. The host processes seen above are different: short-lived clients that MATLAB starts and that end by themselves.
"""
from __future__ import annotations

import os
import shutil
import threading
import time
from pathlib import Path
from typing import NamedTuple

import femm_queue
import paths

from . import FAIL, PASS, Row, fresh, need_file, need_tool, save_log, sha256

TIMEOUT_S = 300
SCRIPT = "check_case"
POLL_S = 0.05


# ----------------------------------------------------------------------------------------------- who is in the job
def _toolhelp():
    """kernel32's process-list and process-time calls, declared (a handle is 64 bits wide), and the PROCESSENTRY32W structure. Built once."""
    global _TOOLHELP
    if _TOOLHELP is None:
        import ctypes
        import types
        from ctypes import wintypes

        class ProcessEntry(ctypes.Structure):
            _fields_ = [("dwSize", wintypes.DWORD), ("cntUsage", wintypes.DWORD), ("th32ProcessID", wintypes.DWORD),
                        ("th32DefaultHeapID", ctypes.c_size_t), ("th32ModuleID", wintypes.DWORD), ("cntThreads", wintypes.DWORD),
                        ("th32ParentProcessID", wintypes.DWORD), ("pcPriClassBase", wintypes.LONG), ("dwFlags", wintypes.DWORD),
                        ("szExeFile", ctypes.c_wchar * 260)]

        k = ctypes.WinDLL("kernel32", use_last_error=True)
        k.CreateToolhelp32Snapshot.argtypes, k.CreateToolhelp32Snapshot.restype = [wintypes.DWORD, wintypes.DWORD], wintypes.HANDLE
        k.Process32FirstW.argtypes = k.Process32NextW.argtypes = [wintypes.HANDLE, ctypes.POINTER(ProcessEntry)]
        k.Process32FirstW.restype = k.Process32NextW.restype = wintypes.BOOL
        k.OpenProcess.argtypes, k.OpenProcess.restype = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD], wintypes.HANDLE
        k.GetProcessTimes.argtypes = [wintypes.HANDLE] + [ctypes.POINTER(wintypes.FILETIME)] * 4
        k.GetProcessTimes.restype = wintypes.BOOL
        k.CloseHandle.argtypes, k.CloseHandle.restype = [wintypes.HANDLE], wintypes.BOOL
        _TOOLHELP = types.SimpleNamespace(ctypes=ctypes, wintypes=wintypes, k=k, ProcessEntry=ProcessEntry)
    return _TOOLHELP


_TOOLHELP = None
PROCESS_QUERY_LIMITED_INFORMATION = 0x1000     # the least access that GetProcessTimes needs: it reads, it cannot end or change a process
FILETIME_AT_UNIX_EPOCH = 116444736000000000    # a FILETIME counts 100 ns from 1601-01-01; this is 1970-01-01


def process_table() -> dict[int, tuple[int, str]]:
    """{pid: (parent pid, program name)} of every process running now: a read-only snapshot (Toolhelp), nothing is opened or touched."""
    t = _toolhelp()
    snapshot = t.k.CreateToolhelp32Snapshot(0x2, 0)                    # TH32CS_SNAPPROCESS
    if snapshot is None or snapshot == t.ctypes.c_void_p(-1).value:
        raise t.ctypes.WinError(t.ctypes.get_last_error())
    table: dict[int, tuple[int, str]] = {}
    try:
        entry = t.ProcessEntry()
        entry.dwSize = t.ctypes.sizeof(entry)
        more = t.k.Process32FirstW(snapshot, t.ctypes.byref(entry))
        while more:
            table[entry.th32ProcessID] = (entry.th32ParentProcessID, entry.szExeFile)
            more = t.k.Process32NextW(snapshot, t.ctypes.byref(entry))
    finally:
        t.k.CloseHandle(snapshot)
    return table


def process_created(pid: int) -> int | None:
    """When process `pid` was created, as a Windows FILETIME, read through a handle opened with PROCESS_QUERY_LIMITED_INFORMATION;
    None if the handle cannot be opened (a protected process, or one that has just ended). With its pid it names a process for
    good: Windows may give the pid of a process that has ended to a new one, but never the same creation time. Read-only."""
    t = _toolhelp()
    handle = t.k.OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION, False, pid)
    if not handle:
        return None
    try:
        created, exited, kernel, user = (t.wintypes.FILETIME() for _ in range(4))
        if not t.k.GetProcessTimes(handle, t.ctypes.byref(created), t.ctypes.byref(exited), t.ctypes.byref(kernel), t.ctypes.byref(user)):
            return None
        return (created.dwHighDateTime << 32) | created.dwLowDateTime
    finally:
        t.k.CloseHandle(handle)


def filetime_now() -> int:
    """The system clock as a FILETIME, the unit of process_created, so that a creation time can be compared with 'now'."""
    return int(time.time() * 10_000_000) + FILETIME_AT_UNIX_EPOCH


class Member(NamedTuple):
    pid: int
    name: str
    first_seen: float
    last_seen: float
    created: int | None = None          # the FILETIME the process was created at; None: its handle could not be opened (known by pid and name)
    ppid: int = 0


class _Record:
    """One member of the job, as the tracker keeps it. `ended_wall` is the clock (FILETIME) of the look at which it was found gone."""
    __slots__ = ("pid", "ppid", "name", "created", "first", "last", "ended_wall")

    def __init__(self, pid: int, ppid: int, name: str, created: int | None, now: float):
        self.pid, self.ppid, self.name, self.created = pid, ppid, name, created
        self.first = self.last = now
        self.ended_wall: int | None = None


class JobTracker:
    """Which processes descend from the one program that run_solver started, and for how long was each seen?

    run_solver's program is a child of this process, started after `baseline` (the process table just before). A process is a
    member if it is such a child, or a child of a member. That is the job's membership as far as it can be read from outside: a
    process started under the job's members belongs to the job, and nothing else can, because the program is created suspended,
    put into the job, and only then run. update(table, now) is called with each new process table; `now` is seconds since the call.

    A process is known by its pid AND its creation time (process_created), not by its pid alone: Windows gives the pid of a
    process that has ended to a new, unrelated process. Three things go wrong with pids alone, and are tested:
      * a new process of the job that happens to get a pid from the baseline was taken for a baseline process and left out;
      * a stranger's child, whose parent pid is the reused pid of a member that has ended, was counted as a member;
      * a stranger's process on the reused pid of a member had the member's "last seen" extended.
    So a member counts as alive only while its pid has ITS creation time; a child of an ended member counts only if it was created
    while that member lived (between the member's creation and the look that found it gone); and a baseline process is only the one
    with the baseline's creation time. If a creation time cannot be read, the process is known by pid and name, and says so.
    `created` is process_created by default; a test passes its own.
    """

    def __init__(self, parent_pid: int, baseline: dict[int, tuple[int, str]], created=process_created):
        self.parent = parent_pid
        self._created = created
        self.before: dict[int, int | None] = {pid: created(pid) for pid in baseline}      # who was there before, and since when
        self._records: list[_Record] = []
        self.updated_at: float | None = None

    # -- who is who
    def _is_baseline(self, pid: int, created: int | None) -> bool:
        if pid not in self.before:
            return False
        old = self.before[pid]
        return old is None or created is None or old == created      # a different creation time: a new process on an old pid

    def was_there_before(self, pid: int) -> bool:
        """True if the process now running as `pid` is one that was already running when the tracker was made."""
        return self._is_baseline(pid, self._created(pid))

    def _alive(self, record: _Record, table: dict[int, tuple[int, str]]) -> bool:
        if record.pid not in table:
            return False
        if record.created is not None:
            now_created = self._created(record.pid)
            if now_created is not None:
                return now_created == record.created                 # the same pid with another creation time is another process
        return table[record.pid][1].lower() == record.name.lower()

    def _born_of_a_member(self, ppid: int, created: int | None) -> bool:
        if ppid == self.parent:
            return True                                                # a new child of the caller: the program run_solver started
        for record in self._records:
            if record.pid != ppid:
                continue
            if record.ended_wall is None:                              # the parent is a member and alive
                if created is None or record.created is None or created >= record.created:
                    return True
            elif created is not None and record.created is not None and record.created <= created <= record.ended_wall:
                return True                                            # an orphan: created while its parent lived
        return False

    def update(self, table: dict[int, tuple[int, str]], now: float, wall: int | None = None) -> None:
        """Take a look: `table` is process_table(), `now` the seconds since the call, `wall` the clock (FILETIME) after the snapshot."""
        wall = filetime_now() if wall is None else wall
        self.updated_at = now
        for record in self._records:                                   # 1. the members that were alive: still the same process?
            if record.ended_wall is None:
                if self._alive(record, table):
                    record.last = now
                else:
                    record.ended_wall = wall
        found = True
        while found:                                                   # 2. new members (a child may be listed before its parent is known: repeat)
            found = False
            parents = {self.parent} | {record.pid for record in self._records}
            for pid, (ppid, name) in table.items():
                if ppid not in parents:
                    continue                                           # the cheap test first: creation times are read only for children of members
                if any(record.pid == pid and record.ended_wall is None for record in self._records):
                    continue                                           # already a member, and alive
                created = self._created(pid)
                if self._is_baseline(pid, created) or any(record.pid == pid and record.created == created and record.ended_wall is not None
                                                           for record in self._records):
                    continue
                if self._born_of_a_member(ppid, created):
                    self._records.append(_Record(pid, ppid, name, created, now))
                    found = True

    def members(self) -> list[Member]:
        return sorted((Member(r.pid, r.name, r.first, r.last, r.created, r.ppid) for r in self._records),
                      key=lambda m: (m.first_seen, m.pid))

    def alive(self, table: dict[int, tuple[int, str]]) -> list[Member]:
        """The members that are running in `table` now (a fresh process_table()), each checked by creation time."""
        return sorted((Member(r.pid, r.name, r.first, r.last, r.created, r.ppid) for r in self._records
                       if r.ended_wall is None and self._alive(r, table)), key=lambda m: (m.first_seen, m.pid))

    def last_seen(self, name: str) -> float | None:
        """When a member of this program name (Windows names are not case-sensitive) was last seen; None if none ever was."""
        seen = [m.last_seen for m in self.members() if m.name.lower() == name.lower()]
        return max(seen) if seen else None

    def still_present(self, now: float) -> list[str]:
        return [m.name for m in self.members() if m.last_seen >= now]

    def summary(self, returned_at: float, program: str = "matlab.exe", label: str = "MATLAB",
                looks: int | None = None, failed: int = 0, last_error=None) -> str:
        """Which processes joined the job, and how long after the last `program` the job was empty (run_solver returned).

        `looks` and `failed` are how many times the process list was taken and how many of those failed: a list that could not be
        read must never read as an empty job."""
        members = self.members()
        warning = ""
        if failed:
            warning = (f"{failed} of {looks if looks is not None else '?'} looks at the process list FAILED"
                       f"{' (' + str(last_error) + ')' if last_error else ''}, so the list below may be incomplete. ")
        if not members:
            if looks == 0 or (failed and failed >= (looks or 0)):
                return warning + "UNKNOWN: no process of the job was seen, and the process list could not be read."
            return warning + f"No process of the job was seen in {looks if looks is not None else 'any'} looks (the program ended before the first look)."
        names = ", ".join(f"{m.name} ({m.first_seen:.1f}-{m.last_seen:.1f} s{'' if m.created is not None else ', creation time unreadable'})" for m in members)
        last = self.last_seen(program)
        text = warning + f"{len(members)} process(es) joined the job: {names}."
        if last is None:
            return text + f" No {program} was seen."
        later = sorted({m.name for m in members if m.last_seen > last})
        text += f" The last {program} was seen at {last:.2f} s; run_solver returned (the job was empty) at {returned_at:.2f} s, {returned_at - last:.2f} s later."
        text += f" Still there after {label}: {', '.join(later)}." if later else f" Nothing else was in the job after {label}."
        return text


class Watcher(threading.Thread):
    """Feeds a JobTracker with the process table every POLL_S seconds until stop() is called. Read-only; it never ends anything.

    It counts its looks and the ones that failed (`looks`, `failed`, `last_error`), so that a process list that cannot be read is
    reported as that and not as a job that held nothing. `error` is set if the tracker itself raised (a bug): the thread then stops."""

    def __init__(self, tracker: JobTracker, zero: float):
        super().__init__(daemon=True, name="job-watcher")
        self.tracker, self.zero = tracker, zero
        self.looks = self.failed = 0
        self.last_error: Exception | None = None
        self.error: Exception | None = None
        self._stop_now = threading.Event()

    def run(self) -> None:
        while not self._stop_now.is_set():
            self.looks += 1
            try:
                self.tracker.update(process_table(), time.perf_counter() - self.zero)
            except OSError as error:                     # the list could not be taken: counted, and said in the summary
                self.failed += 1
                self.last_error = error
            except Exception as error:                   # a bug here must not pass for "nothing was in the job": say so in the summary
                self.error = error
                return
            self._stop_now.wait(POLL_S)

    def stop(self) -> None:
        self._stop_now.set()
        self.join(5)


# ----------------------------------------------------------------------------------------------- the MATLAB script
HEAD = """% check_case.m: written by checks/matlab_check.py of the motor project, run by `matlab -batch check_case`.
% It leaves result.txt (KEY=value lines); an error stops it and makes MATLAB exit with a nonzero code.
fid = fopen('result.txt', 'w');
closer = onCleanup(@() fclose(fid));
fprintf(fid, 'VERSION=%s\\n', version);
fprintf(fid, 'RELEASE=%s\\n', version('-release'));
fprintf(fid, 'MATLABROOT=%s\\n', matlabroot);
try
    fprintf(fid, 'MAXNUMCOMPTHREADS=%d\\n', maxNumCompThreads);
catch
    fprintf(fid, 'MAXNUMCOMPTHREADS=-1\\n');
end
try
"""

PLANT_PART = """    % --- 1. the plant: load a copy, do not simulate it, close it without saving ---
    params = load('motor structure.mat');
    names = fieldnames(params);
    for k = 1:numel(names), assignin('base', names{k}, params.(names{k})); end
    fprintf(fid, 'PARAMS=%d\\n', numel(names));
    assert(~bdIsLoaded('Improved_plant'), 'a model named Improved_plant is already loaded');
    load_system('Improved_plant.slx');
    fprintf(fid, 'PLANT_LOADED=%d\\n', bdIsLoaded('Improved_plant'));
    fprintf(fid, 'PLANT_STOPTIME=%s\\n', get_param('Improved_plant', 'StopTime'));
    close_system('Improved_plant', 0);
    fprintf(fid, 'PLANT_CLOSED=%d\\n', ~bdIsLoaded('Improved_plant'));
"""

TINY_PART = """    % --- 2. a small model that this script builds: Sine -> Gain -> To Workspace ---
    mdl = 'chk_gain'; gain = 3; amp = 2; w = 5; bias = 0.5;
    assert(~bdIsLoaded(mdl), 'the model name chk_gain is already in use');
    new_system(mdl);
    add_block('simulink/Sources/Sine Wave', [mdl '/Sine'], 'Amplitude', num2str(amp), 'Frequency', num2str(w), 'Bias', num2str(bias), 'SampleTime', '0.01');
    add_block('simulink/Math Operations/Gain', [mdl '/Gain'], 'Gain', num2str(gain));
    add_block('simulink/Sinks/To Workspace', [mdl '/OutU'], 'VariableName', 'u_log', 'SaveFormat', 'Array');
    add_block('simulink/Sinks/To Workspace', [mdl '/OutY'], 'VariableName', 'y_log', 'SaveFormat', 'Array');
    add_line(mdl, 'Sine/1', 'Gain/1'); add_line(mdl, 'Sine/1', 'OutU/1'); add_line(mdl, 'Gain/1', 'OutY/1');
    set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'FixedStepDiscrete', 'FixedStep', '0.01', 'StopTime', '1');
    out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
    t = out.tout; u = out.get('u_log'); y = out.get('y_log');
    errGain = max(abs(y - gain * u)); errSine = max(abs(u - (amp * sin(w * t) + bias)));
    fprintf(fid, 'TINY_SAMPLES=%d\\n', numel(t));
    fprintf(fid, 'TINY_ERR_GAIN=%.3e\\n', errGain);
    fprintf(fid, 'TINY_ERR_SINE=%.3e\\n', errSine);
    close_system(mdl, 0);
    fprintf(fid, 'TINY_CLOSED=%d\\n', ~bdIsLoaded(mdl));
    assert(errGain < 1e-12, 'the gain block did not multiply by 3');
    assert(errSine < 1e-12, 'the sine block is not 2 sin(5 t) + 0.5');
"""

TAIL = """    fprintf(fid, 'OK=1\\n');
catch err
    fprintf(fid, 'ERROR=%s | %s\\n', err.identifier, strrep(err.message, newline, ' '));
    rethrow(err);
end
"""


def script_text() -> str:
    return HEAD + PLANT_PART + TINY_PART + TAIL


# ----------------------------------------------------------------------------------------------- the row
def read_result(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if path.is_file():
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            key, _, value = line.partition("=")
            if _:
                values[key.strip()] = value.strip()
    return values


def thread_count(text: str | None, cap: int) -> tuple[int | None, str]:
    """(the number of computational threads MATLAB reported if it is between 1 and `cap`, else None; and the words for the row). `text` is the value of
    MAXNUMCOMPTHREADS in result.txt: None if the script never got that far, -1 if maxNumCompThreads could not be asked."""
    try:
        count = int(text)
    except (TypeError, ValueError):
        return None, f"MATLAB did not report its computational threads (maxNumCompThreads), so the cap of {cap} is not shown to hold"
    if count < 1:
        return None, f"MATLAB did not report its computational threads (maxNumCompThreads gave {count}), so the cap of {cap} is not shown to hold"
    if count > cap:
        return None, f"MATLAB reports {count} computational threads, above the cap of {cap} (at most {cap}: R10)"
    return count, f"MATLAB reports {count} computational thread{'' if count == 1 else 's'} (at most {cap})"


def run(folder: Path):
    matlab = need_tool("matlab")
    plant = need_file("simulink_plant", paths.path("simulink_plant"))
    params = need_file("motor_structure_mat", paths.path("motor_structure_mat"))
    fresh(folder)
    case = "matlab -batch: load a copy of Improved_plant.slx (not simulated), then simulate a small model it builds"
    cap = paths.setting("femm", "workers")
    expected = (f"exit code 0, the script reaches its last line, the references are unchanged, the job empties, MATLAB computes with at most {cap} threads "
                "(R10; -singleCompThread)")
    originals = {path: sha256(path) for path in (plant, params)}
    for path in (plant, params):
        shutil.copy2(path, folder / path.name)
    copies = {folder / path.name: originals[path] for path in (plant, params)}
    (folder / f"{SCRIPT}.m").write_text(script_text(), encoding="utf-8")

    before = process_table()
    tracker = JobTracker(os.getpid(), before)
    zero = time.perf_counter()
    watcher = Watcher(tracker, zero)
    watcher.start()
    try:
        result = femm_queue.run_solver([matlab, "-sd", folder, "-logfile", folder / "matlab.log", "-singleCompThread", "-batch", SCRIPT],
                                       label="check_tools MATLAB -batch", cwd=folder, timeout_s=TIMEOUT_S, threads=4)
        returned_at = time.perf_counter() - zero
    finally:
        watcher.stop()
    seconds = round(result.seconds, 2)                  # from the program's start until its job was empty (a wait for the slot is not in it)
    job = tracker.summary(returned_at, looks=watcher.looks, failed=watcher.failed, last_error=watcher.last_error)
    if watcher.error is not None:
        job += f" (The job watcher stopped early: {type(watcher.error).__name__}: {watcher.error}.)"
    save_log(folder, "matlab", result)
    (folder / "job_members.txt").write_text(job + "\n\n" + "\n".join(
        f"{m.pid:>7}  {m.name:<28} first {m.first_seen:7.1f} s  last {m.last_seen:7.1f} s" for m in tracker.members()) + "\n", encoding="utf-8")

    unchanged = all(sha256(path) == digest for path, digest in originals.items())
    missing_copies = [path.name for path in copies if not path.is_file()]          # a copy that vanished is not "unchanged"
    copies_unchanged = not missing_copies and all(sha256(path) == digest for path, digest in copies.items())
    if not unchanged:
        references = "a reference file CHANGED"
    elif missing_copies:
        references = f"a scratch copy is MISSING ({', '.join(missing_copies)})"
    elif not copies_unchanged:
        references = "a scratch copy CHANGED"
    else:
        references = "references unchanged"
    values = read_result(folder / "result.txt")
    if result.timed_out:
        yield Row("MATLAB", case, expected, f"timed out after {TIMEOUT_S} s and the job was ended. {job}", FAIL, seconds)
        return
    errors = values.get("ERROR") or (f"exit code {result.returncode}" if result.returncode != 0 else "")
    threads, threads_text = thread_count(values.get("MAXNUMCOMPTHREADS"), cap)
    ok = (result.returncode == 0 and values.get("OK") == "1" and unchanged and copies_unchanged and values.get("PLANT_LOADED") == "1"
          and threads is not None)
    measured = (f"exit code {result.returncode}; MATLAB {values.get('RELEASE', '?')}: plant loaded and closed, gain model y = 3u "
                f"(max error {values.get('TINY_ERR_GAIN', '?')}); {threads_text}; {references}"
                + (f"; {errors}" if errors and not ok else "") + f". Job: {job}")
    yield Row("MATLAB", case, expected, measured, PASS if ok else FAIL, seconds)
