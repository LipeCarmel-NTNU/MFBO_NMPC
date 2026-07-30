"""Internals of the MFBO-NMPC pipeline.

You do not run these modules. Run run_pipeline.py at the project root, and edit
run_config.py at the project root to change a setting.

  driver.py            the phase logic: design, proposals, refits, ledger
  matlab_interface.py  the request and result exchange with MATLAB
  provenance.py        manifest, ledger, environment capture, reconciliation
  phi_surrogate.py     the fit of the fidelity surrogate phi(z)
"""
