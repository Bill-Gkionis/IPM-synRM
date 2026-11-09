"""What the femm_session tests share: the magnet they solve, the worker function, and two Windows process look-ups.

The worker function lives here, not in the test file, because a *spawned* worker (the only way Windows starts a
process) begins as a fresh Python that has to import the function it runs. pytest's own test modules are not
importable that way; this plain module is, as long as tests/ is on sys.path (conftest.py arranges that).
"""
import contextlib
import ctypes
import os
from ctypes import wintypes
from pathlib import Path

import femm

import femm_session

# ---------------------------------------------------------------------------------------------- the magnet


def draw_magnet(hc):
    """learn/03's model on the FEMM that open_femm() started: a 10 x 2 mm magnet in a 60 x 60 mm box of air.

    hc is the magnet's coercivity in A/m: 1,030,272 is the N45 of the motor; 0 makes it a plain block of air-like
    material, so nothing drives any flux.
    """
    femm.newdocument(0)                                              # 0 = magnetics problem
    femm.mi_probdef(0, "millimeters", "planar", 1e-8, 20, 30)        # units, 2D planar, 20 mm deep
    femm.mi_getmaterial("Air")
    femm.mi_addmaterial("N45", 1.05, 1.05, hc)                       # name, mu_x, mu_y, coercivity [A/m]
    femm.mi_drawrectangle(-30, -30, 30, 30)                          # the air box
    femm.mi_drawrectangle(-5, -1, 5, 1)                              # the magnet
    for (x, y), material, magdir in (((0, 0), "N45", 90), ((0, 20), "Air", 0)):
        femm.mi_addblocklabel(x, y)
        femm.mi_selectlabel(x, y)
        femm.mi_setblockprop(material, 1, 0, "", magdir, 0, 0)
        femm.mi_clearselected()
    femm.mi_addboundprop("A=0", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # no flux leaves the box
    for x, y in ((0, -30), (30, 0), (0, 30), (-30, 0)):
        femm.mi_selectsegment(x, y)
    femm.mi_setsegmentprop("A=0", 0, 1, 0, 0)
    femm.mi_clearselected()


def solve_magnet(fem_path, hc):
    """Draw the magnet, save it as fem_path, solve it and return By [T] at the magnet's centre.

    Works on the FEMM that femm_session.open_femm() started. Deletes an old .ans first, so a finished .ans
    next to fem_path can only be this run's answer.
    """
    fem_path = Path(fem_path)
    fem_path.parent.mkdir(parents=True, exist_ok=True)
    fem_path.with_suffix(".ans").unlink(missing_ok=True)
    draw_magnet(hc)
    femm.mi_saveas(str(fem_path))
    femm.mi_analyze(1)                                               # 1 = no progress window
    femm.mi_loadsolution()
    return femm.mo_getb(0, 0)[1]


def magnet_worker(fem_path, hc, both_alive, results):
    """One worker process: open its own FEMM, solve the magnet, report (hc, port, FEMM pid, By, worker pid).

    both_alive is a multiprocessing.Barrier(2): each worker waits at it with its FEMM still open, so the two
    FEMMs are certain to exist at the same moment. If this worker fails it breaks the barrier, so the other one
    does not sit waiting for it.
    """
    try:
        session = femm_session.open_femm()
        try:
            by = solve_magnet(fem_path, hc)
            both_alive.wait(timeout=60)
            results.put(("ok", hc, femm_session.port(session), femm_session.pid(session), by, os.getpid()))
        finally:
            femm_session.close_femm()
    except BaseException as exc:
        both_alive.abort()
        results.put(("error", hc, repr(exc)))
        raise


# ------------------------------------------------------------------------------- two Windows look-ups (ctypes)
_TH32CS_SNAPPROCESS = 0x2
_PROCESS_QUERY_LIMITED_INFORMATION = 0x1000


class _ProcessEntry(ctypes.Structure):
    """PROCESSENTRY32W from the Windows API: one row of the list of running processes."""
    _fields_ = [("dwSize", wintypes.DWORD), ("cntUsage", wintypes.DWORD), ("th32ProcessID", wintypes.DWORD),
                ("th32DefaultHeapID", ctypes.c_size_t), ("th32ModuleID", wintypes.DWORD),
                ("cntThreads", wintypes.DWORD), ("th32ParentProcessID", wintypes.DWORD),
                ("pcPriClassBase", wintypes.LONG), ("dwFlags", wintypes.DWORD),
                ("szExeFile", wintypes.WCHAR * 260)]


def _kernel32():
    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    k32.CreateToolhelp32Snapshot.argtypes = [wintypes.DWORD, wintypes.DWORD]
    k32.CreateToolhelp32Snapshot.restype = wintypes.HANDLE
    k32.Process32FirstW.argtypes = [wintypes.HANDLE, ctypes.POINTER(_ProcessEntry)]
    k32.Process32FirstW.restype = wintypes.BOOL
    k32.Process32NextW.argtypes = [wintypes.HANDLE, ctypes.POINTER(_ProcessEntry)]
    k32.Process32NextW.restype = wintypes.BOOL
    k32.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    k32.OpenProcess.restype = wintypes.HANDLE
    k32.GetPriorityClass.argtypes = [wintypes.HANDLE]
    k32.GetPriorityClass.restype = wintypes.DWORD
    k32.GetCurrentProcess.restype = wintypes.HANDLE
    k32.SetPriorityClass.argtypes = [wintypes.HANDLE, wintypes.DWORD]
    k32.SetPriorityClass.restype = wintypes.BOOL
    k32.CloseHandle.argtypes = [wintypes.HANDLE]
    return k32


def process_info(pid):
    """(image name, parent pid) of the running process with this pid, or None if there is no such process."""
    k32 = _kernel32()
    snapshot = k32.CreateToolhelp32Snapshot(_TH32CS_SNAPPROCESS, 0)
    if snapshot is None or snapshot == ctypes.c_void_p(-1).value:
        raise ctypes.WinError(ctypes.get_last_error())
    try:
        entry = _ProcessEntry()
        entry.dwSize = ctypes.sizeof(entry)
        more = k32.Process32FirstW(snapshot, ctypes.byref(entry))
        while more:
            if entry.th32ProcessID == pid:
                return entry.szExeFile, entry.th32ParentProcessID
            more = k32.Process32NextW(snapshot, ctypes.byref(entry))
        return None
    finally:
        k32.CloseHandle(snapshot)


def priority_class(pid):
    """The Windows priority class of a running process, e.g. 0x4000 = BELOW_NORMAL_PRIORITY_CLASS."""
    k32 = _kernel32()
    handle = k32.OpenProcess(_PROCESS_QUERY_LIMITED_INFORMATION, False, pid)
    if not handle:
        raise ctypes.WinError(ctypes.get_last_error())
    try:
        return k32.GetPriorityClass(handle)
    finally:
        k32.CloseHandle(handle)


def set_own_priority_class(priority_class):
    """Set the priority class of this (the test) process, e.g. 0x20 = NORMAL_PRIORITY_CLASS."""
    k32 = _kernel32()
    if not k32.SetPriorityClass(k32.GetCurrentProcess(), priority_class):
        raise ctypes.WinError(ctypes.get_last_error())


@contextlib.contextmanager
def priority_class_restored():
    """Whatever the block does to the priority class of THIS process (femm_session.lower_priority is process-wide: open_femm calls it), the class the process had when
    the block began is put back when it ends. conftest.py wraps every test in this, so a test that opens a FEMM inside the pytest process no longer leaves every test after
    it at below normal (the timing tests among them). Yields that class."""
    before = priority_class(os.getpid())
    try:
        yield before
    finally:
        if priority_class(os.getpid()) != before:
            set_own_priority_class(before)


# ------------------------------------------------------------------------------- the default motor, solved at a few points
def default_model_points(folder, tag, points, fine=False):
    """Open a FEMM of its own, draw the default motor, solve (theta_e, id, iq) points, close the FEMM.

    Returns (list of the solved numbers, SHA-256 of the .ans of the last point). Every name carries `tag`, so workers that run
    this at the same time never share a file. The numbers are the model's own (psi_d, psi_q, T_age, T_wst, ...), untouched.
    A module-level function on purpose: a spawned worker has to import it (see the note at the top of this file).
    """
    import hashlib

    import motorgeo as mg
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    p = mg.Params()
    fem = folder / f"{tag}.fem"
    ans = folder / f"{tag}_keep.ans"
    femm_session.open_femm()
    try:
        mg.to_femm(p, mg.build(p), str(fem), fine=fine)
        femm.mi_close()
        model = mg.FemmModel(p, str(fem))
        rows = []
        for theta, i_d, i_q in points:
            r = model.solve(theta, i_d, i_q, keep_ans=str(ans))
            rows.append({k: r[k] for k in ("theta_e", "id", "iq", "psi_a", "psi_b", "psi_c", "psi_d", "psi_q", "T_age", "T_wst", "T_dq")})
        return rows, hashlib.sha256(ans.read_bytes()).hexdigest()
    finally:
        femm_session.close_femm()


def default_model_worker(args):
    """One worker process: default_model_points() with a tag of its own."""
    folder, tag, points, fine = args
    return default_model_points(folder, tag, points, fine)
