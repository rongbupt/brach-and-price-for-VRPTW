'''
增加时间窗
'''
from math import sqrt
import re

def distance_cal(cor1,cor2):
  return round(((cor1[0]-cor2[0])**2 + (cor1[1]-cor2[1])**2)**0.5,2)

def read_vrplib_file(file_path):
    file = open(file_path, 'r')
    file_name = 'UNknow'
    for line in file:
   
      if line.startswith('C1') or line.startswith('R1') or line.startswith('RC1'):
        file_name = line.split().pop()   
      
      if line.startswith('NUMBER'):
        line = next(file)
        capacity = int(line.split()[1])
        numofvehicles = int(line.split()[0])
      
      if line.startswith("CUST NO"):
        next(file)
        coordinates = {}
        demand = {}
        ready_time = {}
        due_time = {}
        service_time = {}
        line_ = next(file)
        while len(line_.split())==7:
          row = line_.split()
          coordinates[int(row[0])] = (int(row[1]),int(row[2]))
          demand[int(row[0])] = int(row[3])
          ready_time[int(row[0])] = int(row[4])
          due_time[int(row[0])] = int(row[5])
          service_time[int(row[0])] = int(row[6])
          try:
            line_ = next(file)
          except StopIteration:
            break
        break    
    file.close()
    numofnodes = len(coordinates)#包含depot
    dis_mat = []
    for i in range(numofnodes):
      dis_mat_row = []
      for j in range(numofnodes):
        dis_mat_row.append(distance_cal(coordinates[i],coordinates[j]))
      dis_mat.append(dis_mat_row)
    instance_data = (file_name, capacity, numofvehicles, coordinates, demand,\
      ready_time, due_time, service_time,numofnodes,dis_mat)
    return instance_data

# (file_name, capacity, numofvehicles, coordinates, demand,\
#       ready_time, due_time, service_time,numofnodes,dis_mat) = read_vrplib_file("100_customer/c101.txt")
# print("OK")

   



