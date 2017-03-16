class Downscaler(object):
   """
   Downscale large-scale trajectories
   """
   def __init__(self):
      pass

   def generate(self, values):
      """

      Arguments:
         values (np.array): A large-scale trajectory

      Returns:
         np.array: A 4D array (T, X, Y, V)
      """
      raise NotImplementedError()


class Qq(object):
   """
   Quantile mapping between a fine-scale model and a coarse-scale model
   """
   def __init__(self, parameters):
      self.parameters = parameters

   def generate(self, values):
      pass
