class Variable(object):
   def __init__(self, name, units=None, label=None):
      self.name = name
      self.units = units
      if label is None:
         self.label = self.name
      else:
         self.label = label

   def pretty(self):
      if self.units is None:
         label = "%s" % (self.label.capitalize())
      else:
         label = "%s (%s)" % (self.label.capitalize(), self.units)

      label = label.replace("_", " ")
      return label
