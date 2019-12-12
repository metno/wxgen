class Variable(object):
    def __init__(self, name, units="", label=None):
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

    def __hash__(self):
        return hash((self.name, self.units, self.label))

    def __eq__(self, other):
        if other is None:
            return False
        else:
            return self.name == other.name

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return self.name
