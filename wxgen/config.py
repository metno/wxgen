import pandas as pd
import numpy as np

class Config(object):
    """Configuration for constraints (weights) when matching segments"""
    def __init__(self, filename):
        self.points = pd.read_csv(filename, sep=";").sort_values("variable")
        expected_cols = set(['lat', 'lon', 'variable', 'weight'])
        if (actual_cols := set(self.points.columns)) != expected_cols:
            raise NotImplementedError(f"Expected cols {expected_cols} but received {actual_cols}")

    def group_by_variable(self) -> dict[str, pd.DataFrame]:
        """Dict with {variable_name: points_for_this_variable}"""
        return {var: df_points for var, df_points in self.points.groupby("variable")}
    
    def variables(self) -> list[str]:
        return list(self.group_by_variable.keys())
    
    def weights_as_arr(self) -> np.ndarray[float]:
        return self.points["weight"].values