from gurobipy import Model, GRB, quicksum
from collections import defaultdict
from collections import OrderedDict
import numpy as np
import random
import ast
import csv
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def plot_vehicle_stacks(vehicle_stacks, L, W, H):
    # 為每個客戶生成一個隨機顏色
    customer_colors = {}
    for v, stacks in vehicle_stacks.items():
        for stack in stacks:
            customer = stack["客戶"]
            if customer not in customer_colors:
                # 隨機生成 RGB 顏色
                customer_colors[customer] = (random.random(), random.random(), random.random())
    
    for v, stacks in vehicle_stacks.items():
        if not stacks:
            print(f"車輛 {v} 無分配貨物，跳過視覺化。")
            continue
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title(f"車輛 {v} 的貨物堆疊結果")
        
        # 設置車輛容器邊界
        ax.set_xlim([0, L])
        ax.set_ylim([0, W])
        ax.set_zlim([0, H])
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        
        # 繪製每個貨物
        for stack in stacks:
            x, y, z = stack["位置"]
            l, w, h = stack["尺寸"]
            customer = stack["客戶"]
            
            # 獲取客戶的顏色
            color = customer_colors[customer]
            
            # 計算貨物的 8 個頂點
            vertices = [
                [x, y, z],
                [x + l, y, z],
                [x + l, y + w, z],
                [x, y + w, z],
                [x, y, z + h],
                [x + l, y, z + h],
                [x + l, y + w, z + h],
                [x, y + w, z + h]
            ]
            
            # 定義貨物的 6 個面
            faces = [
                [vertices[0], vertices[1], vertices[5], vertices[4]],  # 底面
                [vertices[2], vertices[3], vertices[7], vertices[6]],  # 頂面
                [vertices[0], vertices[3], vertices[7], vertices[4]],  # 左面
                [vertices[1], vertices[2], vertices[6], vertices[5]],  # 右面
                [vertices[0], vertices[1], vertices[2], vertices[3]],  # 前面
                [vertices[4], vertices[5], vertices[6], vertices[7]]   # 後面
            ]
            
            # 繪製貨物
            ax.add_collection3d(Poly3DCollection(faces, facecolors=color, linewidths=1, edgecolors='k', alpha=0.8))
        
        plt.show()

model = Model("ChainOsr")
Z = 4   # 區域數量
N = 12  # 客戶數量
N0 = 49 # 節點數量(包含倉庫和客戶節點)
O = 6   # 旋轉方向
L, W, H = 300,170,165  #貨櫃長寬高
M = 1e6  # 極大值
C0 = 6   #每單位才積委外的費用
vi = {}  #客戶 i 所有貨物的總才積
vit = {} #客戶 i 的第 t 件貨物的才積

lito = {} #客戶 i 的第 t 件貨物以擺放方向 o 時的長
wito = {} #客戶 i 的第 t 件貨物以擺放方向 o 時的寬
hito = {} #客戶 i 的第 t 件貨物以擺放方向 o 時的高
pito = {} #擺放方向指標
fit = {}  #脆弱性貨物指標
Gi = {}   #客戶 i 的貨物數量

Neb = defaultdict(set)  #在非重疊服務區域 b 內的客戶集合
Nobb1 = defaultdict(list) #在服務區域 b 和 b+1 重疊部分內的客戶集合
Nb = defaultdict(list)   #在服務區域b的預定義之車輛路線的節點，沒有倉庫
Ab = defaultdict(list)    #所有已完成路線的弧線集合
base_dir = r"datasets/N12_A4_S20250102" 

with open(os.path.join(base_dir, "customerInfo.csv"), newline='', encoding='utf-8-sig') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        customer_id = int(row['客戶'].strip())
        cargo_count = int(row['貨物數'].strip())
        totalVolume = int(row['總才積'].strip())
        Gi[customer_id] = cargo_count
        vi[customer_id] = totalVolume

with open(os.path.join(base_dir, "goods.csv"), newline='', encoding='utf-8-sig') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        customer_id = int(row['客戶'])
        t = int(row['貨物'])
        cargoVolume = int(row['材積'])
        good_l = int(row['長'])
        good_w = int(row['寬'])
        good_h = int(row['高'])
        fragile = int(row['脆弱性'])  # 0 或 1

        # 定義 6 種方向的 (l, w, h) 排列
        orientation_results = [
            (good_l, good_w, good_h),
            (good_l, good_h, good_w),
            (good_w, good_l, good_h),
            (good_w, good_h, good_l),
            (good_h, good_l, good_w),
            (good_h, good_w, good_l),
        ]

        # 初始化該客戶
        if customer_id not in lito:
            lito[customer_id] = []
            wito[customer_id] = []
            hito[customer_id] = []

        # 拆開 l, w, h 分別加入對應 list
        l_list = [l for l, w, h in orientation_results]
        w_list = [w for l, w, h in orientation_results]
        h_list = [h for l, w, h in orientation_results]

        lito[customer_id].append(l_list)
        wito[customer_id].append(w_list)
        hito[customer_id].append(h_list)

        if customer_id not in pito:
            pito[customer_id] = {}
        pito[customer_id][t] = {}

        for o in range(O):  # 方向1~6
            flag = int(row[f'方向{o+1}'])
            pito[customer_id][t][o] = flag

        if customer_id not in vit:
            vit[customer_id] = {}
        vit[customer_id][t] = cargoVolume
        fit[(customer_id, t)] = fragile

with open(os.path.join(base_dir, "serviceArea.csv"), newline='', encoding='utf-8-sig') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        customer_id = int(row['客戶'])
        service_vector = [
            int(row['服務區域1'].strip()),
            int(row['服務區域2'].strip()),
            int(row['服務區域3'].strip()),
            int(row['服務區域4'].strip()),
        ]
        if service_vector.count(1) == 1:
            area_index = service_vector.index(1) + 1
            Neb[area_index].add(customer_id)
        for b in range(3):  # b = 0, 1, 2 → 對應區域1~3
            if service_vector[b] == 1 and service_vector[b + 1] == 1:
                Nobb1[b + 1].append(customer_id)  

    Neb = OrderedDict((k, sorted(v)) for k, v in sorted(Neb.items()))
    Nobb1 = OrderedDict((k, sorted(v)) for k, v in sorted(Nobb1.items()))
# print(Nobb1)
with open(os.path.join(base_dir, "routes.csv"), newline='', encoding='utf-8-sig') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        area = int(row['區域'])
        for key in row:
            if key != '區域':
                value = row[key].strip()
                if value != '' and value != '0':
                    Nb[area+1].append(int(value))
    Nb = dict(Nb) 
# print(Nb)
with open(os.path.join(base_dir, "routeArcs.csv"), newline='', encoding='utf-8-sig') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        area = int(row['區域'])
        start = int(row['起點'])
        end = int(row['終點'])
        Ab[area+1].append((start, end))  
    Ab = OrderedDict(sorted(Ab.items())) 
# print(Ab)

#定義各區弧線集合
arc_vars = []
for b in Ab:
    for (i, j) in Ab[b]:
        arc_vars.append((i, j, b))
#所有弧線集合
A = set()
for arc_list in Ab.values():
    A.update(arc_list)
# print(A)

#變數定義

zeta= model.addVars(
    list((i, b) for b in Nobb1 for i in Nobb1[b]) +
    list((i, b + 1) for b in Nobb1 for i in Nobb1[b]),
    vtype=GRB.BINARY,
    name="zeta"
)#先定義分配重疊區域內的顧客變數

delta = model.addVars(
    ((i,b) for i in range(1, N + 1) for b in range(1,Z+1)),
    vtype=GRB.BINARY,
    name="delta"
)#顧客被分配到區域變數

psi = model.addVars(
    arc_vars,
    vtype=GRB.BINARY,
    name="psi"
)#arc變數

omega = model.addVars(
    range(1, N + 1),
    vtype=GRB.BINARY,
    name="omega"
)#若客戶節點 i 被分配至外包車

kappa = model.addVars(
    ((i,b) for i in range(1, N + 1) for b in range(1,Z+1)),
    vtype=GRB.BINARY,
    name="kappa"
)#輔助變數


gamma = model.addVars(
    ((i,j) for i in range(1,N + 1) for j in range(1,N + 1) if i != j), 
    vtype=GRB.BINARY, 
    name="gamma") #同車先後順序變數

alpha = model.addVars(
    ((i, t, o) for i in range(1, N + 1) for t in range(1,Gi[i] + 1) for o in range(O)),
    vtype=GRB.BINARY,
    name="rotate_ito"
)#貨物擺放方向

xprime = model.addVars(
    ((i, t, j, tprime) for i in range(1, N + 1) for t in range(1,Gi[i] + 1) for j in range(1, N + 1)  for tprime in range(1,Gi[j] + 1)
    if (i == j and t != tprime) or (i != j)), 
    vtype=GRB.BINARY, 
    name="f"
)#貨物擺在前面

yprime = model.addVars(
    ((i, t, j, tprime) for i in range(1, N + 1) for t in range(1,Gi[i] + 1) for j in range(1, N + 1)  for tprime in range(1,Gi[j] + 1)
    if (i == j and t != tprime) or (i != j)), 
    vtype=GRB.BINARY, 
    name="r" 
)#右邊

zprime = model.addVars(
   ((i, t, j, tprime) for i in range(1, N + 1) for t in range(1,Gi[i] + 1) for j in range(1, N + 1)  for tprime in range(1,Gi[j] + 1)
    if (i == j and t != tprime) or (i != j)), 
    vtype=GRB.BINARY, 
    name="t" 
)# 上面

x = model.addVars(
    ((i, t) for i in range(1, N + 1) for t in range(1,Gi[i] + 1)), 
    vtype=GRB.INTEGER, 
    name="x")  # Box X-coordinate

y = model.addVars(
    ((i, t) for i in range(1, N + 1) for t in range(1,Gi[i] + 1)), 
    vtype=GRB.INTEGER, 
    name="y"
)  # Box Y-coordinate

z = model.addVars(
    ((i, t) for i in range(1, N + 1) for t in range(1,Gi[i] + 1)), 
    vtype=GRB.INTEGER, 
    name="z"
)  # Box Z-coordinate

R = model.addVars(
    range(1, Z + 1),
    vtype=GRB.CONTINUOUS,
    name="R"
) #第b區自有車輛裝載總才積

R_max = model.addVar(vtype=GRB.CONTINUOUS, name="R_max") #所有自有車輛服務區域之最大裝載總才積
R_min = model.addVar(vtype=GRB.CONTINUOUS, name="R_min") #所有自有車輛服務區域之最小裝載總才積

#目標式
obj1 = quicksum(C0 * vi[i] * omega[i] for i in range(1, N+1))   # 第一層目標（較重要）
obj2 = R_max - R_min   

model.ModelSense = GRB.MINIMIZE
model.setObjectiveN(obj1, index=0, priority=2, name="primary_cost")
model.setObjectiveN(obj2, index=1, priority=1, name="balance")

#目標式計算相關
for b in range(1, Z + 1):
    model.addConstr(R_max >= R[b])

for b in range(1, Z + 1):
    model.addConstr(R_min <= R[b])

for b in range(1, Z + 1):
    model.addConstr(
        R[b] == quicksum(
            vi[i] * (delta[i, b] - kappa[i, b])
            for i in Nb[b]
        )
    )

#輔助變數相關
for b in range(1, Z + 1): 
    for i in Nb[b]:
        model.addConstr(kappa[i, b] <= omega[i])
        model.addConstr(kappa[i, b] <= delta[i,b])
        model.addConstr(kappa[i, b] >= omega[i] + delta[i,b] - 1) 
#路線與區域指派相關
#3A
for b in Nobb1:  # Nobb1[b] = list of customers in overlap between b and b+1
    for i in Nobb1[b]:
        model.addConstr(zeta[i, b] + zeta[i, b + 1] == 1)
#3B
for b in Neb:
    for i in Neb[b]:
        model.addConstr(delta[i, b] + omega[i] == 1)
#3C
for b in range(2, Z + 1):  # Z2 ~ Zu
    if (b-1) in Nobb1:
        for i in Nobb1[b-1]:
            model.addConstr(delta[i, b] + omega[i] == zeta[i, b])
            model.addConstr(delta[i, b-1] + omega[i] == zeta[i, b-1])
#3D 不需要
# for i in Nobb1.get(1, []):  # Nobb1[1] 是 1,2 重疊的客戶
#     model.addConstr(delta[i, 1] == zeta[i, 1])

#3E
for b in range(1, Z + 1):
    customers_in_zone = set(Nb[b])  # 該區域的顧客
    for i in range(1, N + 1):
        if i not in customers_in_zone:
            model.addConstr(delta[i, b] == 0)

#3F
for b in Ab:
    for (i, j) in Ab[b]:
        if i != 0 and j != 0:
            model.addConstr(
                psi[i, j, b] <= (delta[i, b] + delta[j, b]) / 2,
            )
#3G
for b in Ab:
    # 建立 incoming mapping：j → 所有 (i, j)
    incoming = defaultdict(list)
    for (i, j) in Ab[b]:
        incoming[j].append((i, j))

    for j in Nb[b]:  # j ∈ N_b
        model.addConstr(
            quicksum(psi[i, j, b] for (i, j) in incoming.get(j, [])) == delta[j, b],
        )

#3H
for b in Ab:
    # 建立 incoming mapping：j → 所有 (i, j)
    incoming = defaultdict(list)
    for (j, k) in Ab[b]:
        incoming[j].append((j, k))

    for j in Nb[b]:  # j ∈ N_b
        model.addConstr(
            quicksum(psi[j, k, b] for (j, k) in incoming.get(j, [])) == delta[j, b],
        )
#3I
for b in Ab:
    if Neb.get(b):  # 非重疊客戶存在
        depot_outgoing = [(i, j) for (i, j) in Ab[b] if i == 0]
        if depot_outgoing:  # 保險起見，檢查是否有 (0,v)
            model.addConstr(
                quicksum(psi[i, j, b] for (i, j) in depot_outgoing) == 1,
            )
#3J
for b in Ab:
    if Neb.get(b):  # 非重疊客戶存在
        depot_outgoing = [(i, j) for (i, j) in Ab[b] if i == 0]
        if depot_outgoing:  # 保險起見，檢查是否有 (0,v)
            model.addConstr(
                quicksum(psi[i, j, b] for (i, j) in depot_outgoing) <= 1,
            )
#3K
for b in Ab:
    from_depot = [(i, j) for (i, j) in Ab[b] if i == 0]
    to_depot = [(i, j) for (i, j) in Ab[b] if j == 0]

    model.addConstr(
        quicksum(psi[i, j, b] for (i, j) in from_depot) ==
        quicksum(psi[i, j, b] for (i, j) in to_depot),
    )

#4
for (i, j) in A:
    if i != 0 and j != 0:
        model.addConstr(
            quicksum(psi[i, j, b] for b in Ab if (i, j) in Ab[b]) <= gamma[i, j]
        )

#gamma設計
for i in range(1, N): 
    for j in range(1, N): 
        for k in range(1, N): 
            if i != j and i != k and j != k: 
                model.addConstr(
                    gamma[i, j] + gamma[j, k] <= 1 + gamma[i, k],
                )

for i in range(1, N): 
    for j in range(1, N):
        if i != j:
           model.addConstr(gamma[i, j] + gamma[j, i] <= 1) 

for i in range(1, N):  
    for j in range(1, N):  
        if i != j:  
            for b in range(1, Z + 1):  
                for bprime in range(1, Z + 1):
                    if b != bprime: 
                        model.addConstr(
                            gamma[i, j] + gamma[j, i] <= 2 - delta[i, b] - delta[j, bprime],
                        )

for i in range(1, N):  
    for j in range(1, N):  
        if i != j:  
            for b in range(1, Z + 1):  
                model.addConstr(
                    gamma[i, j] + gamma[j, i] >= 1 - (2 - delta[i, b] - delta[j, b])
                )
#5A
for i in range(1, N):  
    for j in range(1, N):  
        if i != j:
            for t in range(1, Gi[i] + 1):
                for tprime in range(1, Gi[j] + 1):
                    model.addConstr(
                        xprime[i, t, j, tprime] + xprime[j, tprime, i, t] +
                        yprime[i, t, j, tprime] + yprime[j, tprime, i, t] + 
                        zprime[i, t, j, tprime] + zprime[j, tprime, i, t] <= 3 * (gamma[i, j] + gamma[j, i]),
                    )

#5B
for i in range(1, N):  
    for j in range(1, N):  
        if i != j:
            for t in range(1, Gi[i] + 1):
                for tprime in range(1, Gi[j] + 1):
                    model.addConstr(
                        xprime[j, tprime, i, t] + yprime[i, t, j, tprime] + yprime[j, tprime, i, t] + zprime[j, tprime, i, t] >= gamma[i, j],
                    )

#5C
for i in range(1, N + 1):
    for t in range(1, Gi[i] + 1):
        for tprime in range(1, Gi[i] + 1):
            if t < tprime:
                model.addConstr(
                    xprime[i, t, i, tprime] + xprime[i, tprime, i, t] +
                    yprime[i, t, i, tprime] + yprime[i, tprime, i, t] +
                    zprime[i, t, i, tprime] + zprime[i, tprime, i, t] >= 1,
                )

#5D
for (i, j) in A:
    if i != 0 and j != 0:
        for t in range(1, Gi[i] + 1):
            for tprime in range(1, Gi[j] + 1):
                if fit[i, t] == 1 and fit[j, tprime] == 0:
                    model.addConstr(
                        xprime[i, t, j, tprime] + xprime[j, tprime, i, t] +
                        yprime[i, t, j, tprime] + yprime[j, tprime, i, t] +
                        zprime[i, t, j, tprime] >= gamma[i, j],
                    )

#5E
for i in range(1, N + 1):
    for t in range(1, Gi[i] + 1):
        for tprime in range(1, Gi[i] + 1):
            if t != tprime:
                if fit[i, t] == 1 and fit[i, tprime] == 0:
                    model.addConstr(
                        xprime[i, t, i, tprime] + xprime[i, tprime, i, t] +
                        yprime[i, t, i, tprime] + yprime[i, tprime, i, t] +
                        zprime[i, t, i, tprime] >= 1,
                    )

#6A 擺放方向
for i in range(1, N + 1):
    for t in range(1, Gi[i] + 1):
        model.addConstr(quicksum(alpha[i, t, o] for o in range(O)) == 1)

#6B~D 
for i in range(1, N+1):
    for t in range(1, Gi[i] + 1):
        for j in range(1, N+1):
            for tprime in range(1, Gi[j] + 1):
                if (i == j and t != tprime) or (i != j):
                    model.addConstr(
                        x[i, t] + quicksum(alpha[i, t, o] * lito[i][t-1][o] for o in range(O)) 
                        <= x[j, tprime] + M * (1 - xprime[i, t, j, tprime]),
                    )
                    model.addConstr(
                        y[i, t] + quicksum(alpha[i, t, o] * wito[i][t-1][o] for o in range(O)) 
                        <= y[j, tprime] + M * (1 - yprime[i, t, j, tprime]),
                    )
                    model.addConstr(
                        z[i, t] + quicksum(alpha[i, t, o] * hito[i][t-1][o] for o in range(O)) 
                        <= z[j, tprime] + M * (1 - zprime[i, t, j, tprime]),
                    )
#6E~F 
for i in range(1, N + 1):
    for t in range(1, Gi[i] + 1):
        model.addConstr(
            x[i, t] + quicksum(alpha[i, t, o] * lito[i][t-1][o] for o in range(O)) <= L + M*omega[i]
        )
        model.addConstr(
            y[i, t] + quicksum(alpha[i, t, o] * wito[i][t-1][o] for o in range(O)) <= W 
        )
        model.addConstr(
            z[i, t] + quicksum(alpha[i, t, o] * hito[i][t-1][o] for o in range(O)) <= H 
        )

#26P
for (i, j) in A:
    if i != 0 and j != 0:
        model.addConstr(omega[j] >= omega[i] + gamma[i,j] - 1)

model.optimize()

print("\n=== 各區客戶區域分配結果 ===")
for b in range(1, Z+1):  # 例如 Z = range(1, num_zones+1)
    print(f"區域 {b}:")
    for i in range(1,N+1):  # 所有客戶
        if (i, b) in delta and delta[i, b].X > 0.5:
            print(f"  客戶 {i} 被分配到區域 {b}")

print("\n=== 各區車輛路徑結果 ===")
for b in range(1, Z+1): 
    print(f"區域 {b} 的路徑:")
    for (i, j) in Ab.get(b, []):  # Ab[b] 是該區域的弧線集合
        if (i, j, b) in psi and psi[i, j, b].X > 0.5:
            print(f"  {i} → {j}")

print("\n=== 分區貨物裝載位置結果 ===")
for b in range(1, Z + 1):  # 假設 Z 是區域總數
    print(f"\n【區域 {b}】")
    for i in lito:
        # 確認顧客 i 有被分配到區域 b，且有被裝載（omega[i] == 0）
        if (i, b) in delta and delta[i, b].X > 0.5 and omega[i].X == 0:
            for t in range(len(lito[i])+1):  # 該客戶的貨物數量
                x_val = x[i, t].X if (i, t) in x else None
                y_val = y[i, t].X if (i, t) in y else None
                z_val = z[i, t].X if (i, t) in z else None
                if x_val is not None and y_val is not None and z_val is not None:
                    print(f"客戶 {i} 的貨物 {t}: x={x_val:.1f}, y={y_val:.1f}, z={z_val:.1f}")



# print("\n=== 各區貨物材積總和")
# for b in range(1, Z + 1):
#     if R[b].X > 0.01:  # 避免浮點誤差影響
#         print(f"區域 {b} 的總裝載材積為：{R[b].X:.2f} 單位³")
#     else:
#         print(f"區域 {b} 的總裝載材積為：0 單位³")


print("\n===顧客被分至外包車")
for i in range(1, N + 1):
    if omega[i].X > 0.01:  # 避免浮點誤差影響
        print(f"顧客 {i} 被分配至租用車輛")
    
print("\n===顧客先後順序")
for i in range(1, N + 1):
    for j in range(1, N + 1):
        if i != j:
            if gamma[i,j].X > 0.01:  # 避免浮點誤差影響
                print(f"顧客 {i} {j} 有先後順序")


area_stacks = {b: [] for b in range(1, Z + 1)}
    
for b in range(1, Z + 1):
    for i in range(1, N):  # 客戶節點
        if delta[i, b].x > 0.5:  # 如果客戶分配給車輛 v
            for t in range(1, Gi[i] + 1):  # 客戶 i 的每件貨物
                x_coord = x[i, t].x
                y_coord = y[i, t].x
                z_coord = z[i, t].x
                for o in range(O):
                    if alpha[i, t, o].x > 0.5:  # 確認旋轉方向
                        if o == 0:  # 長在 x，寬在 y，高在 z
                            length = lito[i][t-1][o]
                            width = wito[i][t-1][o]
                            height = hito[i][t-1][o]
                        elif o == 1:  # 寬在 x，長在 y，高在 z
                            length = lito[i][t-1][o]
                            width = wito[i][t-1][o]
                            height = hito[i][t-1][o]
                        elif o == 2:  # 高在 x，長在 y，寬在 z
                            length = lito[i][t-1][o]
                            width = wito[i][t-1][o]
                            height = hito[i][t-1][o]
                        elif o == 3:  # 長在 x，高在 y，寬在 z
                            length = lito[i][t-1][o]
                            width = wito[i][t-1][o]
                            height = hito[i][t-1][o]
                        elif o == 4:  # 寬在 x，高在 y，長在 z
                            length = lito[i][t-1][o]
                            width = wito[i][t-1][o]
                            height = hito[i][t-1][o]
                        elif o == 5:  # 高在 x，寬在 y，長在 z
                            length = lito[i][t-1][o]
                            width = wito[i][t-1][o]
                            height = hito[i][t-1][o]
                        break 
                
                # 記錄貨物位置和尺寸
                area_stacks[b].append({
                    "客戶": i,
                    "貨物": t,
                    "位置": (x_coord, y_coord, z_coord),
                    "尺寸": (length, width, height)
                })
plot_vehicle_stacks(area_stacks, L, W, H)