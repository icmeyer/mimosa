import re
regex = re.compile('fuel.*')
l = ['fuel1','fuel2','mod','fuel203']
matches = [string for string in l if re.match(regex, string)]
print(matches)
print('fuel1' == l)
