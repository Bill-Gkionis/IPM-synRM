"""b5_prepomax_script.py: PrePoMax's half of B5. It runs INSIDE PrePoMax (its embedded Python 3.12, where the object `pmx` is already there), never in P's environment.

What it does, and why it is a separate file
    B5 asks what PrePoMax 3.0.0 can do without a person at its window. Its own command line (PrePoMax.com --help) has a switch for a Python script, `--script`, which is run in a
    session with no window when `-g No` is given; the script drives the program through the API that PrePoMax documents (Help\\Python\\Python_API_Overview.html in its folder: the
    same objects the menus use, and every change is recorded in the model's command history like a click). So the work is split as B6's is: bench/b5_prepomax.py (P's side)
    writes case.json, starts `PrePoMax.com --script b5_prepomax_script.py --scriptArgs case.json -g No -x Yes ...` through bench.run_program, and reads what comes back; this file
    does the PrePoMax part and writes its numbers to answer.json. They share nothing but those two files and PrePoMax's own --report (a machine-readable record of the run).
    Because it runs in another program's Python it imports the standard library only, and it never starts a program.

What it does, in order (each step is recorded in answer.json with whether it worked, so that a failure says which one)
    1. sets the model's unit system (mm, tonne, s) and model space (plane stress): a script run starts with an empty model that has neither, and an import is refused without them.
    2. imports B4's deck (`pmx.import_file`). PrePoMax's importer reads the whole deck: the mesh and its 6-node triangles, the node and element sets, the material, the section, the
       two boundary conditions, the nodal forces (one "Concentrated_force" load for each node), the output requests. It reads two things wrongly or not at all, found on 2026-11-06:
       the *DLOAD CENTRIF line ("Failed to import distributed load": PrePoMax's own export writes the same line, so the line is right, the reader is not) and the thickness of the
       section (the imported section has none, so the solver would use 1 mm, not the stack's 20 mm, and the nodal forces would then be 20 times too strong against the body force).
    3. reads back what the import made (counts of parts, nodes, sets, material, section, boundary conditions, loads, outputs).
    4. replaces the imported section by one of the same material and element set with the thickness (`sections.add_solid(..., thickness=20)`), and adds the centrifugal load
       (`step.loads.add_centrifugal`, about the z axis, the speed in rad/s).
    5. adds an analysis (`pmx.jobs.add`) with its own work folder, at most `cpus` CPUs, the solver from PrePoMax's settings (ccx_static.exe of CalculiX 2.23, which Task 1 set in
       settings.bin); runs it and waits (`job.run(timeout=...)`). The executable the job used is recorded.
    6. opens the result and reads the largest von Mises stress and its node (`FieldExtremes`).
    7. plots the von Mises stress in PrePoMax's own viewport (the window is hidden but the viewport renders) and saves the picture (`pmx.view.save_image`).
    8. writes answer.json, also when a step failed (the exception is in it), and then lets the exception through, so that PrePoMax's exit code says 3 (a script error).

Run (by bench/b5_prepomax.py; by hand for a look)
    PrePoMax.com --script b5_prepomax_script.py --scriptArgs "case.json" -g No -x Yes -w <folder> --report report.json --timeout 600 --onPrompt Fail
Expect (2026-11-06): the importer's two gaps in the log, status OK, the largest von Mises stress of B4's pressure model at 10,000 rpm (33.454 MPa on the 0.1 mm mesh), a few minutes (CalculiX's part is 6 s: answer.json says what each step took).
"""
import json
import os
import sys
import time
import traceback


def names(collection):
    """The names of the items of a PrePoMax collection, as a plain list of text."""
    return [str(getattr(item, "name", item)) for item in collection]


def counted(collection) -> dict:
    """{name: count} of node sets or element sets."""
    return {str(item.name): int(item.count) for item in collection}


def run_step(answer: dict, label: str, function):
    """Run one step of the list, record in answer["steps"] that it worked (or why not) and how long it took, and return its value; an error is recorded and goes on."""
    started = time.perf_counter()
    try:
        value = function()
    except Exception as error:
        answer["steps"].append({"step": label, "ok": False, "seconds": round(time.perf_counter() - started, 2), "error": f"{type(error).__name__}: {error}"})
        raise
    answer["steps"].append({"step": label, "ok": True, "seconds": round(time.perf_counter() - started, 2)})
    return value


def main(case_file: str) -> None:
    with open(case_file, encoding="utf-8") as handle:
        case = json.load(handle)
    answer = {"python": sys.version.split()[0], "case": case_file, "steps": [], "error": None}
    answer_file = case["answer"]

    def step(label, function):
        return run_step(answer, label, function)

    try:
        def properties():
            pmx.model.properties.unit_system = case["unit_system"]
            pmx.model.properties.model_space = case["model_space"]
        step("set the unit system and the model space", properties)
        step("import B4's deck", lambda: pmx.import_file(case["deck"]))
        step_name = case["step"]
        model_step = pmx.model.steps[step_name]

        def read_back():
            parts = [{"name": str(p.name), "type": str(p.part_type), "nodes": int(p.node_count), "elements": int(p.element_count)} for p in pmx.model.parts]
            loads = {}
            for load in model_step.loads:
                loads[str(load.item_type)] = loads.get(str(load.item_type), 0) + 1
            sets = counted(pmx.model.node_sets)
            named = {n: c for n, c in sets.items() if not n.startswith("Auto_created")}       # the importer makes one set of one node for every nodal force
            return {"parts": parts, "node_sets": named, "auto_created_node_sets": len(sets) - len(named), "element_sets": counted(pmx.model.element_sets), "materials": names(pmx.model.materials),
                    "sections": [{"name": str(s.name), "thickness": s.thickness} for s in pmx.model.sections], "steps": names(pmx.model.steps),
                    "boundary_conditions": names(model_step.bcs), "loads": loads, "field_outputs": names(model_step.field_outputs), "history_outputs": names(model_step.history_outputs)}
        answer["imported"] = step("read back what the import made", read_back)

        def section():
            for old in names(pmx.model.sections):
                pmx.model.sections.remove(old)
            pmx.model.sections.add_solid("Section-B4", case["material"], pmx.create.region("ElementSetName", case["element_set"]), thickness=case["thickness"])
            return {"material": case["material"], "element_set": case["element_set"], "thickness": [s.thickness for s in pmx.model.sections]}
        answer["section"] = step("replace the section by one with the thickness", section)

        def centrifugal():
            model_step.loads.add_centrifugal("Spin", pmx.create.region("ElementSetName", case["element_set"]), axis_point=[0, 0, 0], axis_direction=[0, 0, 1],
                                             rotational_speed=case["omega"])
            return {"name": "Spin", "element_set": case["element_set"], "rotational_speed_rad_s": case["omega"]}
        answer["centrifugal"] = step("add the centrifugal load", centrifugal)

        def run_job():
            os.makedirs(case["job_dir"], exist_ok=True)
            job = pmx.jobs.add("Job-1", work_directory=case["job_dir"], number_of_cpus=case["cpus"])
            status = job.run(timeout=case["solver_timeout_s"])
            output = job.get_output() or ""
            return {"name": "Job-1", "status": str(status), "executable": str(job.executable), "argument": str(job.argument), "work_directory": str(job.work_directory),
                    "number_of_cpus": int(job.number_of_cpus), "results_file": str(job.results_file_name), "output_tail": output[-400:]}
        answer["job"] = step("run CalculiX", run_job)
        if answer["job"]["status"] != "OK":
            raise RuntimeError(f"the analysis ended with status {answer['job']['status']}")

        def results():
            pmx.jobs.open_results("Job-1")
            result = pmx.results.current
            last_step = result.steps.last()
            increment = last_step.increments.last()
            peak = increment.fields["STRESS"].get_extremes("MISES")
            return {"result": str(result.name), "steps": len(result.steps), "increments": len(last_step.increments), "fields": [str(n) for n in increment.fields.get_names()],
                    "peak": {"value": float(peak.max_value), "node": int(peak.max_node_id), "minimum": float(peak.min_value), "unit": str(peak.unit)}}
        answer["results"] = step("open the result and read the largest von Mises stress", results)

        def picture():
            pmx.view.current = "results"
            pmx.view.set_result_field("STRESS", "MISES")
            pmx.view.result_type = "ColorContours"
            pmx.view.deformation_type = "TrueScale"             # the true displacement is a few micrometres: the sector looks as it is, and the scale bar is right
            pmx.view.annotations.show_min_max(False, False)     # the arrow from the largest value's label would cross the whole picture: the caption says it
            pmx.view.set_orientation("Front")
            pmx.view.zoom_to_fit()
            written = pmx.view.save_image(case["image"], True, case["magnification"])      # a transparent background: the bench crops what is not the model or the legend
            return {"path": str(written), "bytes": os.path.getsize(written)}
        answer["image"] = step("plot the stress in PrePoMax's own viewport and save the picture", picture)
    except Exception:
        answer["error"] = traceback.format_exc()
        raise
    finally:
        with open(answer_file, "w", encoding="utf-8") as handle:
            json.dump(answer, handle, indent=1)


if "pmx" in globals():                                  # PrePoMax puts `pmx` into the globals of the script it runs; an import (a test) has none and runs nothing
    main(sys.argv[1])
