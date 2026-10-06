"""Generate a simple DIFLOW directional migration map."""

import matplotlib.pyplot as plt
import pandas as pd

from diflow.mapping import plot_directional_map


coordinates = pd.DataFrame(
    {
        "population": ["POP_A", "POP_B", "POP_C"],
        "latitude": [39.0, 39.5, 38.9],
        "longitude": [-84.0, -83.2, -82.7],
    }
)

flows = pd.DataFrame(
    {
        "source": ["POP_A", "POP_B", "POP_B"],
        "destination": ["POP_B", "POP_A", "POP_C"],
        "migration": [0.030, 0.005, 0.018],
        "support": [0.99, 0.71, 0.96],
    }
)

ax = plot_directional_map(
    coordinates,
    flows,
    min_support=0.75,
)
ax.set_title("DIFLOW directional migration map")
plt.tight_layout()
plt.show()
