# Profile-likelihood diagnostics

Bootstrap intervals describe resampling uncertainty, but they do not by
themselves show whether a migration parameter is sharply identifiable.

DIFLOW therefore supports one-dimensional profile likelihoods for either
scaled directional migration parameter.

## Run a profile

After an inference or prepare-only run has saved pairwise spectra:

    diflow profile \
      --spectrum results/spectra/POP_A__POP_B.npy \
      --parameter m_a_to_b \
      --output profile_A_to_B/ \
      --points 15 \
      --starts 10

Or profile the reverse direction with:

    --parameter m_b_to_a

At each fixed migration value, DIFLOW reoptimizes the remaining asymmetric-model
parameters.

Outputs include a CSV profile and PNG/PDF likelihood curves.

## Interpretation

A narrow peak suggests the migration parameter is relatively well identified
under the fitted model.

A broad or flat profile means many migration values produce similar likelihoods,
so the point estimate should be interpreted cautiously.

The reported interval uses the usual one-parameter profile-likelihood chi-square
cutoff as an approximate diagnostic. It is not a substitute for simulation
calibration or bootstrap uncertainty, especially near parameter boundaries.
