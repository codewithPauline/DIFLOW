# Directional migration map

DIFLOW's map layer consumes a standardized edge table rather than raw output
from any one estimator.

## Coordinates

Required columns:

- population
- latitude
- longitude

## Flows

Required columns:

- source
- destination
- migration

Optional columns:

- support
- lower
- upper

Every row is explicitly source -> destination.

## Visual encoding

The first renderer uses:

- arrow direction = migration direction
- arrow width = migration magnitude
- arrow opacity = directional support
- node position = supplied geographic coordinates

Filtering can be applied by minimum migration magnitude or minimum support.

The initial implementation deliberately plots geographic coordinates without
inventing a basemap. A later projected-cartography layer will add real
boundaries, coastlines, watersheds, and publication export while preserving
the same standardized flow table.
