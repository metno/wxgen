class Config(object):
    def __init__(self, filename):
        self.points = list()
        file = open(filename, 'r')
        header = file.readline().strip().split(';')
        Ilat = header.index('lat')
        Ilon = header.index('lon')
        Ivar = header.index('variable')
        if 'weight' in header:
            Iweight = header.index('weight')
        else:
            Iweight = None

        for line in file:
            words = line.strip().split(';')
            curr = dict()
            curr['lat'] = float(words[Ilat])
            curr['lon'] = float(words[Ilon])
            curr['variable'] = words[Ivar]
            if Iweight is not None:
                weight = float(words[Iweight])
            else:
                weight = 1
            curr['weight'] = weight
            self.points += [curr]

        file.close()
