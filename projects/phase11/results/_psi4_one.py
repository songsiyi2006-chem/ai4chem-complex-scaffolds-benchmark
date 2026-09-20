
# One energy per process: this Windows Psi4 build fails PSIO checkpoint
# re-opening when several energies run sequentially inside a single driver.
import sys
import psi4

geom, basis, method, log = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
psi4.set_memory("900 MB")     # stay inside this build's 512 MB - 2 GB band
psi4.set_output_file(log, False)
mol = psi4.geometry("units bohr\n" + geom)
psi4.set_options({"reference": "rhf", "scf__fail_on_maxiter": False})
print("ENERGY", psi4.energy(method + "/" + basis, molecule=mol))
