import importlib.util, json, sys
import numpy as np
spec = importlib.util.spec_from_file_location("eng", "run_phase13_metalloenzyme_pcet_engine.py")
eng = importlib.util.module_from_spec(spec); spec.loader.exec_module(eng)
eng.TIER = "production"
eng.PSI4_PY = eng.discover_psi4_python()
eng.WORKER_PATH.write_text(eng.PSI4_WORKER_SRC)

res = {}
for m in (3, 5):
    r = eng.run_psi4_job(f"a1_ferryl_m{m}", eng.GEOM_FERRYL, 1, m,
                         [["B3LYP", "def2-SVP"]], 2700,
                         density_grid={"box": 4.2, "step": 0.30})
    res[m] = r

e3 = res[3].get("E_accepted"); e5 = res[5].get("E_accepted")
print("B3LYP S=1:", e3, "| B3LYP S=2:", e5)
if e3 and e5 and res[3].get("converged") and res[5].get("converged") and abs(e5 - e3) < 1.5:
    master = json.load(open(eng.RESULTS / "phase13_results.json"))
    a = master["module_13A"]
    a["dE_hilo_eV"] = float((e3 - e5) * eng.HARTREE_EV)
    a["dE_hilo_tier"] = ["B3LYP", "def2-SVP"]
    a["yamaguchi_J_eV"] = float((e3 - e5) / 4.0)
    a["ladder"]["5"] = {**res[5], "S2_sz_estimate": a["ladder"]["5"].get("S2_sz_estimate")}
    a["ladder"]["3"] = {**res[3], "S2_sz_estimate": a["ladder"]["3"].get("S2_sz_estimate")}
    a["dE_hilo_note"] = "same-tier B3LYP pair from the targeted finalize pass"
    json.dump(master, open(eng.RESULTS / "phase13_results.json", "w"),
              indent=1, ensure_ascii=False, default=float)
    print("PATCHED dE_hilo_eV =", a["dE_hilo_eV"])
else:
    print("pair not landed; master record left as-is")
