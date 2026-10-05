`obs.txt`, `lens.txt`, and `probs.txt` are a compact inference capture from the
pinned `hmmlearn==0.3.3` public `GaussianHMM.predict_proba` API. `model.json`
serializes the model parameters from MACS3 3.0.5's upstream
`test/test_HMMR_HMM.py` training fixture. The model has deliberately asymmetric
transition probabilities and full, correlated covariance matrices. The ten
observations are grouped into independent sequences of lengths 3, 4, and 3.

The capture was generated in the pinned MACS3 oracle environment with NumPy
2.5.3, SciPy 1.18.1, scikit-learn 1.9.1, and hmmlearn 0.3.3. Cargo tests read
these recorded files and do not call Python at runtime.
