"""Stand-in solver programs for the femm_queue tests. Not a test, and not a real solver: a few lines that behave like one in the one
way a test looks at (they start children, announce themselves, sleep, or exit at once). They need nothing but the standard library.

    python stand_in_solver.py <mode> <folder> [seconds]

Modes (the ones a test starts; the others are what these start in turn):
    tree      the program starts a child; the program and the child announce themselves (program.pid, child.pid) and sleep
    orphan    the program starts a child; the child starts a grandchild and EXITS, so the grandchild is an orphan (its parent has
              gone, which is the case taskkill /T cannot reach); program and grandchild announce themselves and sleep
    launcher  starts a child that announces itself (child.pid), sleeps <seconds> and writes child.done, and EXITS AT ONCE. The child
              holds none of the launcher's pipes, so nothing but the job tells run_solver that it is still there. This is what
              matlab.exe does when it starts MATLAB.exe
    sleeper   announces itself (program.pid) and sleeps for <seconds> (default 120)

A process announces itself by writing "<its pid> <its creation time>" into <folder>/<name>.pid. The tests open a handle on it at
that moment (it is certainly running, and it is theirs), check the creation time through the handle against the announced one, and
end it only through that handle. A stand-in lives 120 s at most (far longer than any test
waits), so even one that no handle ever reached ends by itself soon.
"""
import os
import pathlib
import subprocess
import sys
import time


def creation_time() -> int:
    """When this process was created, as a Windows FILETIME. With its pid it identifies this process for good, even after the pid is
    given to somebody else: a test that opens a handle checks the time through the handle against this one."""
    import ctypes
    from ctypes import wintypes
    k = ctypes.WinDLL("kernel32", use_last_error=True)
    k.GetCurrentProcess.restype = wintypes.HANDLE
    k.GetProcessTimes.argtypes = [wintypes.HANDLE] + [ctypes.POINTER(wintypes.FILETIME)] * 4
    k.GetProcessTimes.restype = wintypes.BOOL
    created, exited, kernel, user = (wintypes.FILETIME() for _ in range(4))
    if not k.GetProcessTimes(k.GetCurrentProcess(), ctypes.byref(created), ctypes.byref(exited), ctypes.byref(kernel), ctypes.byref(user)):
        raise ctypes.WinError(ctypes.get_last_error())
    return (created.dwHighDateTime << 32) | created.dwLowDateTime


def announce(folder: str, name: str) -> None:
    pathlib.Path(folder, f"{name}.pid").write_text(f"{os.getpid()} {creation_time()}")


def start(*arguments: str, quiet: bool = False) -> None:
    """Start this same script again with other arguments, as a child that is not waited for."""
    options = dict(stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) if quiet else {}
    subprocess.Popen([sys.executable, __file__, *arguments], **options)


def main(mode: str, folder: str, seconds: float) -> None:
    if mode == "tree":
        start("leaf", folder, "child")
        announce(folder, "program")
        print("the program started", flush=True)                # what a solver prints: run_solver has to keep it if it ends the run
        time.sleep(seconds)
    elif mode == "orphan":
        start("middle", folder)
        announce(folder, "program")
        print("the program started", flush=True)
        time.sleep(seconds)
    elif mode == "launcher":
        start("leaf-then-done", folder, str(seconds), quiet=True)
        print("the launcher exits now", flush=True)
    elif mode == "sleeper":
        announce(folder, "program")
        time.sleep(seconds)
    elif mode == "leaf":                                        # started by "tree": [leaf, folder, name]
        announce(folder, sys.argv[3])
        time.sleep(120)
    elif mode == "middle":                                      # started by "orphan": starts the grandchild and ends at once
        start("leaf", folder, "grandchild")
    elif mode == "leaf-then-done":                              # started by "launcher": sleeps, then proves it ran to its end
        announce(folder, "child")
        time.sleep(seconds)
        pathlib.Path(folder, "child.done").write_text("done")
    else:
        raise SystemExit(f"unknown mode {mode!r}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[1] != "leaf" else 120.0)
