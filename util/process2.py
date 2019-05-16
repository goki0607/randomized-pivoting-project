import glob, os
import csv
import math

files = []
for file in glob.glob("*.csv"):
    files = files + [file]

header = ["Pivot Choice","Average CPU Time","# of Times Optimal Solution Found","# of Times Max Iterations Reached","# of Times Lowest Objective Value Reached","# of Times Failed"]
pivots = ["Bland", "Danztig", "Steepest Edge", "Random Edge", "Random Facet", "Clarkson"]
#b = [0,0,0,0,0]
#d = [0,0,0,0,0]
#s = [0,0,0,0,0]
#e = [0,0,0,0,0]
#f = [0,0,0,0,0]
#c = [0,0,0,0,0]

w = 5
h = 6
b = [[0 for x in range(w)] for y in range(h)] 
x = 0
for file in files:
  
  with open(file, mode='r') as cvs_file:
    print(file)
    rows = csv.reader(cvs_file, delimiter=',')
    i = 0
    m = math.inf
    for row in rows:
      if i != 0:
        if row[1] != "-Inf" or row[1] != "Inf":
          if float(row[1]) < m:
            m = float(row[1])
      i += 1
    i = 0
    csv_file = open(file, mode='r')
    rows = csv.reader(csv_file, delimiter=',')
    for row in rows:
      if i != 0:
        i -= 1
        if int(row[4]) == 1:
          b[i][4] += 1
        else:
          b[i][0] += float(row[0])
          if int(row[3]) == 0:
            b[i][1] += 1
          if int(row[3]) == -1:
            b[i][2] += 1
          if float(row[1]) == m:
            b[i][3] += 1
        i += 1
      i += 1
  x += 1
print(x)

i = 0
j = 0
for row in b:
  for col in row:
    b[i][j] = str(b[i][j])
    j += 1
  i += 1
  j = 0

print(b)

new_file = "summary.csv"
write_file = open(new_file, mode='w')
csv_writer = csv.writer(write_file, delimiter=',')
csv_writer.writerow(header)
i = 0
for row in b:
  csv_writer.writerow([pivots[i]]+row)
  i += 1
