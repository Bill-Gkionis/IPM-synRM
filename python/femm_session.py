"""femm_session: one private FEMM 5.0 per Python process.

How a script uses it (every script in this folder that needs FEMM does exactly this):

    import femm
    import femm_session

    femm_session.open_femm()            # a hidden FEMM on its own port, at below-normal priority
    try:
        femm.newdocument(0)             # the ordinary femm.mi_*() / femm.mo_*() calls now reach THAT FEMM
        ...
    finally:
        femm_session.close_femm()       # quits it now. Always: see "A FEMM goes with its script" below

It replaces the pair the scripts used with FEMM 4.2: `femm.openfemm(1)` followed by `mg.lower_femm_priority()`.
open_femm() also returns the session, so `s = femm_session.open_femm()` allows `s.mi_addnode(0, 0)` calls too.

Why not just call femm.openfemm(1), which FEMM 5.0's own binding offers? Three reasons.

1. The port.
   FEMM 5.0 is a small server. It listens on a TCP port, and Python drives it by sending one line of text per
   command. Think of the port as a telephone number. openfemm() always dials 48015, which is also the number an
   interactive FEMM window answers on, and if somebody already answers there it does not start anything: it joins the
   call. So a script run while an interactive FEMM window is open would draw into that window, and the closefemm() at its
   end sends quit() to it, unsaved work and all. A FEMM also answers only one caller at a time (the next one is
   turned away with "femm-qt is already connected to another client"), and a headless FEMM quits when its caller
   hangs up. So eight workers that each call openfemm() at once would not each get a FEMM: one would, and the
   other seven would be turned away.
   open_femm() asks the operating system for a free port (free_port()), starts a new FEMM there, and checks that
   the FEMM it ended up with is one it started itself (pid() is not None). It does not connect to anything else,
   short of an unrelated program taking the number inside FEMM's 0.25 s start-up. That ends in a loud error or a
   retry (the _proc check below refuses any session whose FEMM we did not start, so nothing it answers is ever used),
   with one exception: a stranger that accepts the connection and then never answers would hang the binding's very
   first call, because callfemm waits for a reply with no timeout of its own. That case is unlikely, and open_femm()
   does not close it. The net that catches a hang is femm_queue.femm_pool: its start timeout ends a worker that is
   stuck in open_femm(), and its task_timeout_s ends one that is stuck in a task, each with its whole process tree.
   A script that calls open_femm() with no pool has no such net. Separately, open_femm() even tests a port by trying to
   bind it, not by connecting, because a connection that hangs up again would quit an idle FEMM. (To bind is to claim a
   port number for a socket of its own: the operating system refuses if someone already holds the number, and that refusal
   is how a port can be tested without talking to its owner. _port_taken() says what the test can and cannot see.)
   close_femm() quits only the FEMM that open_femm() started.

2. The priority, set before FEMM starts.
   A priority class is Windows' ranking of a process when several want the CPU at once (normal, below normal,
   idle, ...). A new process does NOT simply get its creator's class: Windows hands the class on only from a
   BELOW_NORMAL or IDLE creator, and otherwise the child starts at NORMAL. (Checked on this PC, 2026-11-05, with a
   throw-away parent of each class: ABOVE_NORMAL gave a NORMAL child, NORMAL a NORMAL child, BELOW_NORMAL a child at
   0x4000.) The class is also handed on only at the moment of creation: lowering the parent later changes nothing
   for children that already exist. That is exactly why lowering first works: lower_priority() puts this process
   at BELOW_NORMAL, so the FEMM started after it is born at BELOW_NORMAL too, which a NORMAL parent could not have
   given it (test_femm_runs_below_normal in tests/test_femm_session.py shows both halves). Asking for a class at the
   moment the FEMM is created is not open to us, because pyfemm-qt starts it with a plain subprocess.Popen, so
   inheriting is the way. FEMM 4.2 could not work this way at all: COM started femm.exe somewhere else, so the
   scripts had to hunt for it afterwards (a PowerShell call, about 0.5 s). The whole Python process is lowered, on
   purpose: it mostly waits for FEMM, and what is to be kept fast is the desktop.

3. The process tree belongs to the script, so there is nothing to hunt for.
   FEMM 4.2 was started by COM, outside the script's process tree, which is why the scripts had to find
   femm.exe and fkn.exe by name afterwards. FEMM 5.0 is an ordinary child: Python starts femm_headlessw.exe (the
   build without a console window), and that starts a short-lived fkn.exe for every solve. Both inherit the class
   from the one lowering before the launch (checked: fkn.exe runs at 0x4000, below normal, too). One FEMM is
   one femm_headlessw.exe: one pid to look at (pid()), to count (learn/04_workers.py counts them), and to end.

A FEMM goes with its script. A headless FEMM quits when the connection to it drops, so it also goes when its Python
process ends, even if that process was killed (checked). close_femm() in a finally still matters: a pool worker lives
on through many tasks, and the next open_femm() in the same process needs the old FEMM closed first.

This is the only module that touches the binding's private parts: femm._default (the session that the module-level
femm.mi_*() functions talk to), and a session's _proc (the Popen of the FEMM it launched), _port and _sock.

Where FEMM is. femm_dir() takes the folder that holds femm_headlessw.exe from config.toml first: the line `femm` in its
[tools] section, read through paths.py, which is the one file that names every program of these scripts. Only if the config
names no FEMM does femm_dir() use FEMM_DIR, the setting the FEMM 5.0 installer wrote. A shell that was opened before
FEMM 5.0 was installed does not have FEMM_DIR (environment variables are copied into a process when it starts), so that
fall-back also reads it from the registry. A config line that names a folder that is not FEMM's is an error, not a reason
to try FEMM_DIR: a stale line should be found and corrected, not hidden behind another FEMM. One more thing follows from
the order: the binding (femm.FemmSession) looks in the FEMM_DIR of the environment first, so open_femm() puts the folder
femm_dir() chose into it, and the FEMM that starts is the one femm_dir() named, whatever the shell had in FEMM_DIR.
"""
from __future__ import annotations

import os
import socket
import subprocess
import sys
import threading
from pathlib import Path

import femm

import paths                          # config.toml: [tools] femm is FEMM's folder

DEFAULT_PORT = 48015                  # where FEMM listens by default, and where femm.openfemm() always goes
_FEMM_PROGRAMS = ("femm_headlessw.exe", "femm_headless.exe", "femm_headless")   # any one of them marks FEMM's folder
_BELOW_NORMAL = 0x4000                # Windows BELOW_NORMAL_PRIORITY_CLASS
_IDLE = 0x40                          # Windows IDLE_PRIORITY_CLASS (lower still)

_lock = threading.Lock()              # open_femm() and close_femm() change state that the whole process shares
_session = None                       # the session open_femm() made in this process ...
_owner = None                         # ... and the pid of the process that made it


# ----------------------------------------------------------------------------- where FEMM is
def _femm_dir_setting() -> tuple[str | None, str | None]:
    """The first FEMM_DIR found, and where: this process's environment, then (Windows) the registry."""
    value = os.environ.get("FEMM_DIR")
    if value:
        return value, "this process's environment"
    if sys.platform == "win32":
        import winreg
        places = ((winreg.HKEY_CURRENT_USER, r"Environment", r"HKCU\Environment"),
                  (winreg.HKEY_LOCAL_MACHINE, r"SYSTEM\CurrentControlSet\Control\Session Manager\Environment",
                   r"HKLM\SYSTEM\CurrentControlSet\Control\Session Manager\Environment"))
        for root, subkey, label in places:
            try:
                with winreg.OpenKey(root, subkey) as key:
                    value, kind = winreg.QueryValueEx(key, "FEMM_DIR")
            except OSError:                                     # no such key or value
                continue
            if kind == winreg.REG_EXPAND_SZ:                    # a value like %ProgramFiles%\femm\bin
                value = winreg.ExpandEnvironmentStrings(value)
            if isinstance(value, str) and value:
                return value, label
    return None, None


def _config_names_femm() -> bool:
    """True if config.toml has a [tools] femm line that is not empty (an empty line, no line, no [tools] section: False)."""
    try:
        return bool(paths.setting("tools", "femm"))
    except KeyError:
        return False


def _holds_femm(folder: Path) -> bool:
    """True if the folder holds one of the programs that mark FEMM's folder (_FEMM_PROGRAMS)."""
    return any((folder / name).is_file() for name in _FEMM_PROGRAMS)


def femm_dir() -> Path:
    """The folder that holds FEMM 5.0's programs: config.toml's [tools] femm first, then FEMM_DIR.

    1. config.toml (paths.py), the line `femm` of [tools]. If it is not empty it is the answer, whatever FEMM_DIR says, and
       it has to be right: a folder that is not there, or that holds none of femm_headlessw.exe, femm_headless.exe and
       femm_headless, raises FileNotFoundError that names the config. It is not passed over for FEMM_DIR, because that
       would hide a stale line behind some other FEMM.
    2. FEMM_DIR, only if the config names no FEMM (an empty line, no line, no [tools] section). Looked up in the
       environment of this process; if the shell predates FEMM 5.0 it is missing there, so the registry is read next
       (on Windows): HKCU\\Environment, then HKLM\\SYSTEM\\CurrentControlSet\\Control\\Session Manager\\Environment.
       The folder must hold femm_headlessw.exe (or femm_headless on other systems).
    Raises FileNotFoundError, naming config.toml's [tools] femm or else FEMM_DIR, if FEMM is not found.
    """
    if _config_names_femm():
        folder = paths.tool("femm")                     # FileNotFoundError, naming the tool, if the folder is not there
        if not _holds_femm(folder):
            raise FileNotFoundError(
                f"[tools] femm in {paths.CONFIG} is {str(folder)!r}, but that folder holds none of "
                f"{', '.join(_FEMM_PROGRAMS)}. Correct that line to the folder that holds femm_headlessw.exe "
                "(for FEMM 5.0 that is its bin folder).")
        return folder
    value, where = _femm_dir_setting()
    if value is None:
        raise FileNotFoundError(
            "FEMM_DIR is not set: not in this process's environment, and not (on Windows) in the registry "
            "under HKCU\\Environment or HKLM\\SYSTEM\\CurrentControlSet\\Control\\Session Manager\\Environment, "
            "and config.toml names no FEMM ([tools] femm). Set femm in config.toml, or FEMM_DIR, to the folder "
            "that holds femm_headlessw.exe (for FEMM 5.0 that is its bin folder).")
    folder = Path(value)
    if not _holds_femm(folder):
        raise FileNotFoundError(
            f"FEMM_DIR is {value!r} (from {where}), but that folder holds none of {', '.join(_FEMM_PROGRAMS)}. "
            "Set FEMM_DIR to the folder that holds femm_headlessw.exe.")
    return folder


def femm_program() -> Path:
    """The file of the FEMM program that open_femm() starts (the first of _FEMM_PROGRAMS that femm_dir() holds)."""
    folder = femm_dir()
    return next(folder / name for name in _FEMM_PROGRAMS if (folder / name).is_file())


_version_dll = None                   # version.dll with its three calls declared; made once by _version_calls()


def _version_calls():
    """version.dll, with the argument and return types of the three calls file_version() makes declared, once.

    Why they must be declared. Unless it is told otherwise, ctypes treats every argument it has to convert as a C `int`
    (32 bits) and every result as a C `int`. A pointer on 64-bit Windows is 64 bits wide. So an address passed as a plain
    number does not fit (ctypes refuses it), and, worse, an address handed back as a result is silently cut to its lower
    32 bits: the program then reads from an address that is not the one the call meant, which crashes or reads the wrong
    bytes without any error. `argtypes` and `restype` say "this one is a pointer" (c_void_p), "this is a DWORD", and
    ctypes converts correctly. They are properties of the function, not of one call, so they are set here once and kept.
    """
    global _version_dll
    if _version_dll is None:
        import ctypes
        from ctypes import wintypes
        dll = ctypes.WinDLL("version", use_last_error=True)
        dll.GetFileVersionInfoSizeW.argtypes = [wintypes.LPCWSTR, ctypes.POINTER(wintypes.DWORD)]
        dll.GetFileVersionInfoSizeW.restype = wintypes.DWORD
        dll.GetFileVersionInfoW.argtypes = [wintypes.LPCWSTR, wintypes.DWORD, wintypes.DWORD, ctypes.c_void_p]
        dll.GetFileVersionInfoW.restype = wintypes.BOOL
        dll.VerQueryValueW.argtypes = [ctypes.c_void_p, wintypes.LPCWSTR, ctypes.POINTER(ctypes.c_void_p),
                                       ctypes.POINTER(wintypes.UINT)]
        dll.VerQueryValueW.restype = wintypes.BOOL
        _version_dll = dll
    return _version_dll


def file_version(path: Path) -> str | None:
    """The version a Windows program file states about itself, for example "5.0.20261101"; None if it states none.

    Every .exe can carry a *version resource*: a block of text inside the file (company, product, file version) that
    Explorer shows under Properties > Details. Reading it takes three calls into version.dll: how big the block is,
    the block itself, and then the one string we want. The string's path inside the block depends on a language and
    code-page pair that the block lists first, so that pair is read first and used to build the path. This is what
    PowerShell's (Get-Item file).VersionInfo.FileVersion shows too, without starting PowerShell.
    Returns None on other systems than Windows, for a file with no version resource and for a file that is not there.
    """
    if sys.platform != "win32":
        return None
    import ctypes
    import struct
    from ctypes import wintypes
    dll = _version_calls()
    size = dll.GetFileVersionInfoSizeW(str(path), None)          # 0 = no such file, or no version resource
    if not size:
        return None
    block = ctypes.create_string_buffer(size)
    if not dll.GetFileVersionInfoW(str(path), 0, size, block):
        return None
    found, length = ctypes.c_void_p(), wintypes.UINT()
    if not dll.VerQueryValueW(block, "\\VarFileInfo\\Translation", ctypes.byref(found), ctypes.byref(length)) or length.value < 4:
        return None
    language, codepage = struct.unpack("<HH", ctypes.string_at(found.value, 4))
    if not dll.VerQueryValueW(block, f"\\StringFileInfo\\{language:04x}{codepage:04x}\\FileVersion",
                              ctypes.byref(found), ctypes.byref(length)) or not found.value:
        return None
    return ctypes.wstring_at(found.value, length.value).rstrip("\0").strip() or None


def femm_version() -> str | None:
    """The version of the FEMM that open_femm() would start, from its program file: "5.0.20261101" for FEMM 5.0.

    This is the installed FEMM, not a statement about files that were made earlier (the flux page's builder compares
    the date of its data with the program's before it names a version). Raises FileNotFoundError, like femm_dir(),
    if FEMM is not found (neither config.toml nor FEMM_DIR leads to it).
    """
    return file_version(femm_program())


# ----------------------------------------------------------------------------- port and priority
def free_port() -> int:
    """A TCP port nobody is using right now: bind to port 0 and let the OS pick, then give the port back.

    To bind is to claim a port number for a socket of its own (module docstring, section 1). Asking for port 0 means
    "any free number": the operating system picks one, getsockname() tells us which, and closing the socket gives
    it back. There is a small gap between giving it back and FEMM taking it. open_femm() copes with a loss there by
    trying once more on a new port.
    """
    while True:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.bind(("127.0.0.1", 0))
            number = s.getsockname()[1]
        if number != DEFAULT_PORT:            # Windows picks above 49151, Linux can pick 48015: never that one
            return number


def _port_taken(number: int) -> bool:
    """True if something already holds this port. Tries to bind it, which tells us without talking to the owner.

    To bind is to claim a port number for a socket of its own (module docstring, section 1). The operating system
    refuses (OSError) if someone holds the number already, so a refusal means "taken" and nobody had to be called.
    Never test a port by connecting to it: a FEMM that has no client yet takes the first connection as its client,
    and FEMM 5.0 quits when its client hangs up. A probe that connects and leaves would quit a FEMM that is not ours.

    Which owners this probe sees. It sees a listener on 127.0.0.1 (exclusive or not), and a listener on 0.0.0.0 (every
    address of the PC) that holds the port exclusively. FEMM 5.0, a Qt server, is the first kind: netstat shows it
    at 127.0.0.1:<port>, so a FEMM that is already running is seen. It does not see a plain wildcard listener, one
    that listens on 0.0.0.0 without claiming the port exclusively (on Windows a second socket may still bind
    127.0.0.1 beside it), nor an IPv6 one (::1 or ::): the probe asks about IPv4 only. All of this was checked on
    this PC (Windows 11, 2026-11-05) with a throw-away listener of each kind, and with a FEMM.

    Why 127.0.0.1: it is where femm.py connects (FemmSession's host), so it is the address that has to be free.

    The second net. If an owner the probe cannot see answers on the number anyway, FemmSession connects to it
    instead of starting a FEMM, so no process of ours stands behind the session (_proc is None). _launch() notices,
    hangs up and raises, and open_femm() tries a new port: a loud error or a retry, not a script that goes on talking
    to the wrong program. (By then FemmSession has sent that owner its one start-up command, setcurrentdirectory,
    and waited for the answer. That and the hang-up are all it does. An owner that accepts the connection and never
    answers would hang that wait, which has no timeout: see section 1 of the module docstring for the net that catches
    it.)
    """
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        try:
            s.bind(("127.0.0.1", number))
        except OSError:
            return True
    return False


def lower_priority() -> None:
    """Put this whole process at below-normal priority. Processes it starts afterwards inherit that.

    Windows: SetPriorityClass(this process, BELOW_NORMAL_PRIORITY_CLASS = 0x4000) through ctypes; a process that
    is already idle stays idle (this never raises a priority). Elsewhere: os.nice() up to niceness 5.
    """
    if sys.platform == "win32":
        import ctypes
        from ctypes import wintypes
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetCurrentProcess.restype = wintypes.HANDLE
        kernel32.GetPriorityClass.argtypes = [wintypes.HANDLE]
        kernel32.GetPriorityClass.restype = wintypes.DWORD
        kernel32.SetPriorityClass.argtypes = [wintypes.HANDLE, wintypes.DWORD]
        kernel32.SetPriorityClass.restype = wintypes.BOOL
        me = kernel32.GetCurrentProcess()               # a pseudo-handle meaning "this process"
        if kernel32.GetPriorityClass(me) == _IDLE:
            return
        if not kernel32.SetPriorityClass(me, _BELOW_NORMAL):
            raise ctypes.WinError(ctypes.get_last_error())
    else:
        niceness = os.getpriority(os.PRIO_PROCESS, 0)
        if niceness < 5:
            os.nice(5 - niceness)


# ----------------------------------------------------------------------------- start and stop FEMM
def _mine():
    """The session this very process opened, else None. (A forked child inherits the variable, not the FEMM.)"""
    return _session if _owner == os.getpid() else None


def _reap(proc: subprocess.Popen, grace: float) -> None:
    """Wait up to `grace` seconds for a process that we started to exit; kill it if it will not."""
    try:
        proc.wait(timeout=grace)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait(timeout=5)


def _hang_up(session: femm.FemmSession) -> None:
    """Drop our connection WITHOUT sending quit() (session.close() would send it)."""
    sock, session._sock = getattr(session, "_sock", None), None
    if sock is not None:
        sock.close()


def _discard(session: femm.FemmSession) -> None:
    """Clean up after a start that failed half way: stop the FEMM process that we launched, and only that one."""
    proc = getattr(session, "_proc", None)
    if proc is None:                                    # nothing of ours was started: never send quit() to a stranger
        _hang_up(session)
        return
    connected = getattr(session, "_sock", None) is not None
    try:
        session.close()                                 # sends quit() if connected, then closes the socket
    except Exception:
        pass
    _reap(proc, grace=2 if connected else 0)


def _launch(number: int) -> femm.FemmSession:
    """Start one FEMM on this port and connect to it. Raises, with nothing of ours left running, if that fails."""
    if _port_taken(number):
        raise RuntimeError(f"port {number} is already in use, by something that is not a FEMM we started")
    # Build the object in two steps. If starting FEMM fails half way, __init__ raises and the object would be lost
    # while the FEMM process it launched keeps running; holding the object lets us stop that process.
    session = femm.FemmSession.__new__(femm.FemmSession)
    try:
        femm.FemmSession.__init__(session, port=number, headless=True)
    except BaseException:
        _discard(session)
        raise
    # The second net, behind _port_taken(): FemmSession tries to connect before it starts anything, and when that works
    # it never starts a FEMM, so _proc stays None. That means somebody else answers on this number: either it took the
    # number after our probe, or it is an owner that _port_taken() cannot see. Hang up (without quit()) and refuse;
    # open_femm() then tries another port.
    if session._proc is None:
        _hang_up(session)
        raise RuntimeError(f"joined a FEMM that was already running on port {number}")
    return session


def open_femm(*, below_normal: bool = True, port: int | None = None) -> femm.FemmSession:
    """Start this process's own FEMM and make it the target of the module-level femm.mi_*() / mo_*() calls.

    Order matters: FEMM_DIR is set to femm_dir()'s answer (the binding looks there first, so it has to name the FEMM
    that config.toml chose, even if the shell's own FEMM_DIR names another folder), then (if below_normal) this
    process is lowered with lower_priority(), and only then is FEMM started, so FEMM inherits the low priority.
    FEMM listens on `port`, or on a free one (free_port()) if none is given; if the start fails, or something
    answers there already, it tries once more on a new free port, so `port` is a request and port(session) says
    which port was obtained. A FEMM that was already running is not joined, sent a command or closed, short of an unrelated
    program taking the number inside FEMM's 0.25 s start-up; that ends in a loud error or a retry, except for a stranger
    that accepts the connection and never answers, which hangs this call (the binding waits for a reply with no
    timeout). Inside femm_queue.femm_pool a hang is ended by the pool's start timeout; called on its own, this has no net.
    (_port_taken() and the _proc check in _launch() are the two nets; their docstrings say what each sees.)

    Returns the session. Raises RuntimeError if this process already has a FEMM open: the module-level calls can
    reach only one, and threads share a process, so give each worker its own process (learn/04_workers.py).
    """
    global _session, _owner
    with _lock:
        if _mine() is not None:
            raise RuntimeError(
                "this process already has a FEMM open: close_femm() it first. The module-level femm.mi_*() calls "
                "reach one FEMM per process, so a second would leave the first one running. Threads share a "
                "process; give every worker its own process (see learn/04_workers.py).")
        folder = femm_dir()                         # a clear error now, before anything is started
        os.environ["FEMM_DIR"] = str(folder)        # femm.FemmSession looks for the programs here FIRST: it must name femm_dir()'s folder
        if below_normal:
            lower_priority()                        # BEFORE the launch: FEMM inherits it
        for attempt in (1, 2):
            number = port if (attempt == 1 and port is not None) else free_port()
            try:
                session = _launch(number)
                break
            except Exception:
                if attempt == 2:
                    raise
        femm._default = session                     # the module-level femm.mi_*() functions talk to this one
        _session, _owner = session, os.getpid()
        return session


def close_femm() -> None:
    """Quit this process's FEMM, wait for it to go (killing it after 5 s if it will not), and forget it.

    Quits only the session that open_femm() opened in this process. Safe to call twice, or with none open.
    """
    global _session, _owner
    with _lock:
        session = _mine()
        _session = _owner = None
        if session is None:
            return
        if femm._default is session:
            femm._default = None
        proc = session._proc
        session.close()                             # sends quit() and closes the socket
        if proc is not None:
            _reap(proc, grace=5)


def pid(session: femm.FemmSession) -> int | None:
    """The process id of the FEMM behind a session. None if it joined a FEMM that was already running (open_femm()
    refuses such a session and never returns one), or once the session has been closed."""
    proc = session._proc
    return None if proc is None else proc.pid


def port(session: femm.FemmSession) -> int:
    """The TCP port the session talks to."""
    return session._port
