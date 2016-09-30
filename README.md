# Large-scale weather generator

This python package generates arbitrarily long time-series of weather data, by sampling from
a database of shorter time-series.

## Database
The database contains  N trajectories  with T days in each trajectory and V number of atmospheric variables.

## Method
Longer trajectories are created by concatinating the shorter trajectories from the database. This is done by matching the end state of one trajectory with the beginning state of another. The matching is done using a specified metric, such as the sum of the square differences  between states (with some kind of normalization strategy  as each atmospheric variable has different variances).
