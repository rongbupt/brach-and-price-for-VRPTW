import pandas as pd
import numpy as np
from gurobipy import *
from collections import namedtuple
import matplotlib.pyplot as plt
from collections import defaultdict
import time
import random
import copy
import matplotlib.pyplot as plt
import itertools
from random import sample
from readdata import read_vrplib_file
from collections import defaultdict

(file_name, C, numofvehicles, coordinates, demand,ready_time, due_time, service_time,numnodes,dis_mat) \
= read_vrplib_file("in/c101.txt")
distance = copy.deepcopy(dis_mat)
distance.append(copy.deepcopy(distance[0]))
for i in distance:
    i.append(copy.deepcopy(i[0]))
coordinates[numnodes] = coordinates[0]
demand[numnodes] = demand[0]
ready_time[numnodes] = ready_time[0]
due_time[numnodes] = due_time[0]
service_time[numnodes] = service_time[0]


#生成单一节点初始路径
def generate_initial_solution():
    path_list = []
    for i in range(1,numnodes):
        path_list.append([i])        
    return path_list

routes_repair = generate_initial_solution()

def get_route_length(route,numnodes):#返回route的长度和该条路径对每个点的访问次数
   #numnodes包含了起点 
    ksit  = [0] * (numnodes)#numnodes包含0点，但是不包含终点
    lastcustomer = 0
    route_length = 0
    for i,cust in enumerate(route):
        ksit[cust] +=1
        route_length += distance[cust][lastcustomer]
        lastcustomer = cust
    route_length += distance[0][lastcustomer]      
    return route_length,ksit

def get_route_length_ori(route,numnodes):
   #numnodes包含了起点 
    ksit  = [0] * (numnodes)#numnodes包含0点
    lastcustomer = 0
    route_length = 0
    for i,cust in enumerate(route):
        ksit[cust] +=1
        route_length += dis_mat[cust][lastcustomer]
        lastcustomer = cust
    route_length += dis_mat[0][lastcustomer]      
    return route_length,ksit

def get_pi():
    pi =[]
    pi = [c.Pi for c in  master.getConstrs()]
    pi.insert(0, 0)#加入起点
    pi.append(0)#加入终点
    return pi

def get_reduced_cost_mat():#包含了起点和终点
    reduced_cost_mat = []
    for i in range(len(pi)):
        reduced_cost_mat_temp = []
        for j in range(len(pi)):
            if(i==j):
                reduced_cost_mat_temp.append(0)
            else:
                reduced_cost_mat_temp.append(distance[i][j]-pi[i]/2-pi[j]/2)
        reduced_cost_mat.append(reduced_cost_mat_temp)
#     reduced_cost_mat.append(reduced_cost_mat[0])
#     for i in reduced_cost_mat:
#         i.append(i[0])
    return reduced_cost_mat

            
def get_neighbor():
    neighbor = []
    for i in range(numnodes):
        neighbor_row = [numnodes]
        for j in range(numnodes):
            if(ready_time[i]+service_time[i]+distance[i][j] <= due_time[j] and pi[j]>0 and i!=j):
                neighbor_row.append(j)
        neighbor.append(neighbor_row)
    neighbor.append([])
    return neighbor

def get_sort_reduced_cost_mat():
    sort_reduced_cost_mat = []
    for i in range(len(pi)):
        sort_reduced_cost_mat.append(sorted(range(len(pi)), key=lambda k: -reduced_cost_mat[i][k]))
#         sort_reduced_cost_mat.append(sorted(range(1,len(pi)), key=lambda k: -distance[i][k]))
    return sort_reduced_cost_mat
def get_neighbor_2():
    #减半稀疏
    neighbor = []
    for i in range(len(pi)):
        neighbor_temp = []
        for j in sort_reduced_cost_mat[i][-param_collection:]:
            if(j!=i):
                neighbor_temp.append(j)
        # neighbor_temp.insert(random.randint(0,5),1)
#         neighbor_temp.append(1)
        neighbor_temp.insert(0,Terminate)
        neighbor.append(neighbor_temp)
    return neighbor

def C_carry(path):
    result = 0
    for i in path:
        result += demand[i]
    return result

def re_cost(path):
    front = 0
    L = 0
    for i in range(1,len(path)):
        L += reduced_cost_mat[front][path[i]]
        front = path[i]
    return L

def pulse(vi,rp,qp,tp,p,root,t_root):
    if isFeasible(vi,qp,tp):
        if checkBounds(vi,tp,rp,root,t_root,p):
            if rollback(vi,tp,rp,p):
                p_new = copy.deepcopy(p)
                p_new.append(vi)
                qp_new = qp+demand[vi]
                for vj in neighbor[vi]:
                    if(vj not in p_new):
                        rp_new = rp+reduced_cost_mat[vi][vj]
                        tp_new = max(ready_time[vj],tp+service_time[vi]+distance[vi][vj])
                        pulse(vj,rp_new,qp_new,tp_new,p_new,root,t_root)
    #print(root,"  ",t_root,"已最优")

def isFeasible(vi,qp,tp):
    if(demand[vi]+qp>C or tp>due_time[vi]):
        #print('超时')
        return False
    return True

def checkBounds(vi,tp,rp,root,t_root,p):#root表示当前计算bounding的起始点,t_root表示当前计算bonding时到达root的时间
    global path_result
    if(vi==Terminate and rp>primal_bound[root,t_root]):
        #print("到终点，不更新primal_bound") 
        return False
    if(vi==Terminate and rp<0):
        if(root==0 and t_root==0):
            path_result_list[root,t_root].append(copy.deepcopy(p))
    if(vi==Terminate and rp<=primal_bound[root,t_root]):
        primal_bound[root,t_root] = rp
        calculated[root,t_root] = True
        path_result[root,t_root] = copy.deepcopy(p)
        if(root==0 and t_root==0):
            path_result_list[root,t_root].append(copy.deepcopy(p))
#         if(len(path_result_list[0,0])>=1):
#             return False            
        #print("到终点，更新primal_bound ",primal_bound[root,t_root])
        return False
    time_step_left = due_time[0] - (int((due_time[0]-tp)/deta)+1)*deta
    time_step_right = due_time[0] - int((due_time[0]-tp)/deta)*deta
    #print(time_step_left,time_step_right)
    if ((vi,time_step_right) in calculated and rp+primal_bound[vi,time_step_right]<primal_bound[root,t_root]):
        p_temp = copy.deepcopy(p)
        p_temp.extend(path_result[vi,time_step_right])
        if(C_carry(p_temp)<=C): 
            primal_bound[root,t_root] = rp+primal_bound[vi,time_step_right]
            calculated[root,t_root] = True
            path_result[root,t_root] = copy.deepcopy(p_temp)
            if(root==0 and t_root==0):
                path_result_list[root,t_root].append(copy.deepcopy(p_temp))
    if(time_step_left<ready_time[0] or (vi,time_step_left) not in calculated) and (time_step_right not in calculated):
        return True
    if((vi,time_step_left) in calculated) and ((vi,time_step_right) in calculated):
        bound_right = primal_bound[vi,time_step_right]+(time_step_right-tp)*naive_bound
        bound_left = primal_bound[vi,time_step_left]        
        bound_value = max(bound_left,bound_right)
        if(rp+bound_value > primal_bound[root,t_root]):
        #print(vi,"超出bounds")
            return False
        else:
            return True
    if((vi,time_step_left) in calculated) and ((vi,time_step_right) not in calculated):
        if(rp+primal_bound[vi,time_step_left] > primal_bound[root,t_root]):
        #print(vi,"超出bounds")
            return False
        else:
            return True
    if(time_step_left<ready_time[0] or (vi,time_step_left) not in calculated) and ((vi,time_step_right) in calculated):
        bound_right = primal_bound[vi,time_step_right]+(time_step_right-tp)*naive_bound      
        if(rp+bound_right > primal_bound[root,t_root]):
        #print(vi,"超出bounds")
            return False
        else:
            return True


def most_eff():
    most_efficient = 100000000
    for i in range(len(reduced_cost_mat)-1):
        for j in neighbor[i]:
            temp = max(ready_time[j]-due_time[i],service_time[i]+distance[i][j])
            if(temp > 0):
                if reduced_cost_mat[i][j]/temp < most_efficient:
                    most_efficient = reduced_cost_mat[i][j]/temp
    return most_efficient



def rollback(vi,tp,rp,p):
    if(len(p)<=1):
        return True
    else:
        prevNode = p[-1]
        directNode = p[-2]
        directCost = rp-reduced_cost_mat[prevNode][vi]-reduced_cost_mat[directNode][prevNode]+reduced_cost_mat[directNode][vi]
        if(directCost <= rp ):
            return False
    return True

#deta:step_size
def bounding():
    T = copy.deepcopy(due_time[0])
    T -= deta
    while T>ready_time[0]:
        for vi in range(numnodes):
            p = []
            rp = 0
            qp = 0
            tp = max(ready_time[vi],T)
            root = vi
            t_root = T
            pulse(vi,rp,qp,tp,p,root,t_root)
        T -= deta

#设置某条边是forbid
def edge_forbid(forbid_from,forbid_to,distance):
    distance[forbid_from][forbid_to] = large_number

#设置某条边是set
def edge_set(set_from,set_to,distance):
    for i in range(1,len(distance)):
        if(i!=set_to):
            distance[set_from][i] = large_number
    # for i in range(1,len(distance)):
    #     if(i!=set_from):
    #         distance[i][set_to] = large_number

#从route中提取边
def get_edge_from_route(route):
    result = []
    for i in range(len(route)-1):
        result.append((route[i],route[i+1]))
    return result

#判断一条路径是否应该被ban
def route_ban_ornot(route,forbid_from,forbid_to):
    # if((forbid_from,forbid_to) in get_edge_from_route(route)):# or ((forbid_to,forbid_from) in get_edge_from_route(route)):
    #     return True
    # else:
    #     return False
    if(forbid_from in route or forbid_to in route):
        return True
    else:
        return False
def route_set_ornot(route,set_from,set_to):
#     for i in range(len(route)-1):
#         if(route[i]==set_from and route[i+1]!=set_to) or (route[i]==set_to and route[i+1]!=set_from):
#             return True
#     return False
    # if (set_from in route and set_to not in route) or (set_to in route and set_from not in route):
    #     return True
    # else:
    #     return False
    if(set_from in route or set_to in route):
        return True
    else:
        return False

#选择要分支的边
def choose_edge_to_branch(routes,coeff):
    max_temp = 100
    index_chosen = 100
    for i in range(len(routes)):
        if(max_temp > (coeff[i]-0.6)**2):
            max_temp = (coeff[i]-0.6)**2
            index_chosen = i
    return get_edge_from_route(routes[i])[0]

#设计分支节点
class branching_node:
    def __init__(self,node_distance,node_routes_pool,index,father,depth):
        self.int_solution = -1 #当前的一个整数解，用来改善上界
        self.continues_solution = -1 #当前的连续解，也是分支下界
        self.node_distance = node_distance #当前的距离矩阵
        self.node_routes_pool = node_routes_pool #当前的routes_pool
        self.node_y_continues = []#当前连续解的y值，用以确定分支
        self.node_routes_continues = []#当前连续解的routes，用以确定分支
        self.node_routes_int = []#原始整数解
        self.index = index
        self.father = father
        self.depth = depth
        self.node_lower_bound = 0
        self.child_0 = -1
        self.child_1 = -1
        self.forbidden = []
        self.setted = []
    def show(self):
        print('int_solution: ')
        print(self.int_solution)
        print('continues_solution: ')
        print(self.continues_solution)
        print('node_lower_bound: ')
        print(self.node_lower_bound)
        print('node_y_continues: ')
        print(self.node_y_continues)
        print('node_routes_continues: ')
        print(self.node_routes_continues)
        print('node_routes_int: ')
        print(self.node_routes_int)

#更新节点的属性
def update_attribute(node):
    #input self.node_distance; self.node_routes_pool
    #更新 self.int_solution; self.continues_solution; self.node_y_continues; self.node_routes_continues 
    global master
    global routes_pool
    global distance
    global pi
    global reduced_cost_mat
    global sort_reduced_cost_mat
    global neighbor
    global deta
    global Terminate
    global primal_bound
    global calculated
    global T
    global T_list
    global naive_bound
    global path_result
    global path_result_list
    distance = node.node_distance
    routes_pool = node.node_routes_pool
    initialroutecount = len(routes_pool)
    ksitranspose = defaultdict(list)
    route_length    = defaultdict(float)
    #把初始解转换成RMP需要的变量形式
    for r in range(initialroutecount):
        route_length[r],ksitranspose[r] = get_route_length(routes_pool[r],numnodes)

    master = Model("MASTER")
    master.modelSense = GRB.MINIMIZE

    y = {}
    for r in range(initialroutecount):
        y[r] = master.addVar(lb=0.0, vtype=GRB.CONTINUOUS,obj=route_length[r],name='y_%s' % (r))
    master.update()
    custconstr = {}
    for i in range(1,numnodes):#约束中从1号点开始，不包含起点和终点
        custconstr[i] = master.addConstr(
          quicksum(ksitranspose[r][i] * y[r] for r in range(initialroutecount)) >= 1,
                   'cust_%s' % (i))
    master.Params.OutputFlag = 0
    master.update()
    master.optimize()

    param_iteration = 200
    iteration_cout = 0
    start_time = time.clock()
    # while iteration_cout<param_iteration:
    column_generation_result = []
    while True:
        pi = get_pi()
        reduced_cost_mat = get_reduced_cost_mat()
        neighbor = get_neighbor()
        path_result = {}
        path_result_list = defaultdict(list)
        deta = int(due_time[0]/10)
        Terminate = 101
        primal_bound = {}
        calculated = {}
        T = copy.deepcopy(due_time[0])
        T_list = []
        naive_bound = most_eff()
        while T>ready_time[0]:
            T_list.append(T)
            for vi in range(numnodes):
                primal_bound[vi,T] = 10000000000
            primal_bound[Terminate,T] = 0
            T -= deta
        T = copy.deepcopy(due_time[0])
        primal_bound[0,0] = 0
        bounding()
        pulse(0,0,0,0,[],0,0)
        result_temp = path_result_list[0,0]

        new_path = [i[1:] for i in result_temp]
        if(primal_bound[0,0]>=-0.1):
            break
        # if(len(column_generation_result)>15 and column_generation_result[-1]==column_generation_result[-10]):
        #     break

        K = len(routes_pool)
        ksitranspose = defaultdict(list)
        route_length = defaultdict(float)

        for i,j in enumerate(new_path):
            route_length[K],ksitranspose[K] = get_route_length(j,numnodes)
            routes_pool.append(j)

            col = Column()
            for i in range(1,numnodes):
                col.addTerms(ksitranspose[K][i], custconstr[i])
            y[K] = master.addVar(lb=0.0, vtype=GRB.CONTINUOUS,obj=route_length[K], column=col,name='y_%s' % (K))
            master.update()     
            K +=1

        master.Params.OutputFlag = 0
        master.update()
        master.optimize()
        print(iteration_cout,"  ",master.Objval)
        column_generation_result.append(master.Objval)
        iteration_cout += 1
    end_time = time.clock()
    print("当前点列生成时长----",end_time-start_time)

    node.continues_solution = master.ObjVal
    node.node_lower_bound = master.ObjVal
    basic_index = []
    y_continues = []
    for i in y:
        if(y[i].x!=0):
            basic_index.append(i)
            y_continues.append(y[i].x)
    node.node_y_continues = y_continues
    node.node_routes_continues = copy.deepcopy([i for i in np.array(routes_pool)[basic_index]])
    NumberofVariable = len(routes_pool)
    for i in range(NumberofVariable):
        y[i].vType = GRB.BINARY
    master.update()   
    master.optimize()
    node.int_solution = master.ObjVal
    print("整数解:    ",master.ObjVal)
    path_index = []
    for i in y:
        if(y[i].x !=0):
            path_index.append(i)
    node.node_routes_int = copy.deepcopy([i for i in np.array(routes_pool)[path_index]])
    

#抬全局下界
def global_lowerbound_update(index):
    while index>1:
        child_0 = branching_node_dic[branching_node_dic[index].father].child_0
        child_1 = branching_node_dic[branching_node_dic[index].father].child_1
        branching_node_dic[branching_node_dic[index].father].node_lower_bound = min(branching_node_dic[child_0].node_lower_bound,\
                                                                                    branching_node_dic[child_1].node_lower_bound)
        index = branching_node_dic[index].father


r1 = generate_initial_solution()
branching_node_list = []
branching_node_dic = {}
branching_node_dic[0] = branching_node([],[],0,-1,0)
branching_node_list.append(branching_node_dic[0])
branching_node_dic[1] = branching_node(distance,r1,1,-1,1)
branching_node_list.append(branching_node_dic[1])
global_upperbound = 100000000
global_lowerbound = 0
node_index = 1
large_number = 100000000
global_upperbound_routes = []

def branching_procedure():
    global node_index
    global global_upperbound
    global global_lowerbound
    global large_number
    global global_upperbound_routes
    start_time = time.clock()
    while len(branching_node_list)>1:
        current_node = branching_node_list.pop()
        print("processing "+"node "+str(current_node.index))
        update_attribute(current_node)
        if(current_node.int_solution < global_upperbound):
            global_upperbound = current_node.int_solution
            global_upperbound_routes = current_node.node_routes_int
        if(current_node.index % 2 == 0):
            global_lowerbound_update(current_node.index)
        if(current_node.int_solution - current_node.continues_solution < 0.5 or current_node.continues_solution >= global_upperbound):
            print("********CUTTING*******")
            continue
        if(current_node.depth == param_depth_limit):
            print("*****depth_limit*****")
            continue
        (branching_from,branching_to) = choose_edge_to_branch(current_node.node_routes_continues,current_node.node_y_continues)
        
        node_index += 1
        distance_new = copy.deepcopy(current_node.node_distance)
        routes_pool_new = copy.deepcopy([i for i in current_node.node_routes_continues])
        routes_pool_new.extend(copy.deepcopy([i for i in current_node.node_routes_int]))
        edge_forbid(branching_from,branching_to,distance_new)
        #edge_forbid(branching_to,branching_from,distance_new)
        for i in routes_pool_new[::-1]:
            if route_ban_ornot(i,branching_from,branching_to):
                routes_pool_new.remove(i)
        routes_pool_new.extend(copy.deepcopy(routes_repair))
        branching_node_dic[node_index] = branching_node(distance_new,routes_pool_new,node_index,\
                                                        current_node.index,current_node.depth+1)
        branching_node_dic[node_index].forbidden = copy.deepcopy(current_node.forbidden)
        branching_node_dic[node_index].setted = copy.deepcopy(current_node.setted)
        branching_node_dic[node_index].forbidden.append((branching_from,branching_to))
        branching_node_list.append(branching_node_dic[node_index])
        current_node.child_0 = node_index
        
        node_index += 1
        distance_new = copy.deepcopy(current_node.node_distance)
        routes_pool_new = copy.deepcopy([i for i in current_node.node_routes_continues])
        routes_pool_new.extend(copy.deepcopy([i for i in current_node.node_routes_int]))
        edge_set(branching_from,branching_to,distance_new)
        # # edge_set(sett[1],sett[0])#不能反着也加一遍
        for i in routes_pool_new[::-1]:
            if route_set_ornot(i,branching_from,branching_to):
                routes_pool_new.remove(i)
        routes_pool_new.extend(copy.deepcopy(routes_repair))
        branching_node_dic[node_index] = branching_node(distance_new,routes_pool_new,node_index,\
                                                        current_node.index,current_node.depth+1)
        branching_node_dic[node_index].forbidden = copy.deepcopy(current_node.forbidden)
        branching_node_dic[node_index].setted = copy.deepcopy(current_node.setted)
        branching_node_dic[node_index].setted.append((branching_from,branching_to))
        branching_node_list.append(branching_node_dic[node_index])
        current_node.child_1 = node_index
    end_time = time.clock()
    print("全部运行时间----",end_time-start_time)

param_depth_limit = 4#树深度
branching_procedure()

print('global_upperbound: ')
print(global_upperbound)
print('________________________')

print('global_upperbound_routes: ')
print(global_upperbound_routes)
print('________________________')

print("下界：")
print(branching_node_dic[1].node_lower_bound)
print('________________________')
print("根节点连续解")
print(branching_node_dic[1].continues_solution)
print('________________________')


####测试一组解
path_test = global_upperbound_routes
path_draw = []
for i in path_test:
    routes_temp = copy.deepcopy(i)
    routes_temp.insert(0,0)
    routes_temp.append(101)
    path_draw.append(routes_temp)

for i in path_draw:
    plt.plot(np.array([i for i in coordinates.values()])[:,0][i],np.array([i for i in coordinates.values()])[:,1][i],'b-')

plt.plot(np.array([i for i in coordinates.values()])[:,0],np.array([i for i in coordinates.values()])[:,1],'r.')
for i in coordinates:
    plt.text(coordinates[i][0],coordinates[i][1],i,size=12)
plt.show()