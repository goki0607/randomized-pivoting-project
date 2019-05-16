import glob, os
import csv
files = []
for file in glob.glob("*.csv"):
    files = files + [file]

header = ["Pivot Choice","CPU Time","Function Value", "Iterations", "Exit Flag", "Failure"]
pivots = ["Bland", "Danztig", "Steepest Edge", "Random Edge", "Random Facet", "Clarkson"]

for file in files:
  with open(file, mode='r') as cvs_file:
    print(file)
    new_file = file[0:-4] + "2.csv"
    write_file = open(new_file, mode='w')
    csv_writer = csv.writer(write_file, delimiter=',')
    csv_reader = csv.reader(cvs_file, delimiter=',')
    i = -1
    for row in csv_reader:
      if i == -1:
        csv_writer.writerow(header)
      else:
        csv_writer.writerow([pivots[i]]+row)
      i += 1