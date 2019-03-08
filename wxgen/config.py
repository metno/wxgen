class Config(object):
   def __init__(self, filename):
      self.points = list()
      file = open(filename, 'r')
      header = file.readline().strip().split(';')

      for line in file:
         words = line.strip().split(';')
         curr = dict()
         curr['lat'] = float(words[header.index('lat')])
         curr['lon'] = float(words[header.index('lon')])
         curr['variable'] = float(words[header.index('variable')])
         if 'weight' in header:
            weight = float(words[header.index('weight')])
         else:
            weight = 1
         curr['weight'] = weight
         self.points += [curr]

      file.close()
