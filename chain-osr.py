from gurobipy import Model, GRB, quicksum
from collections import defaultdict
from collections import OrderedDict
import numpy as np
import random
import ast
import csv
import math
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
import os

def read_generated_coords(filepath):
    coords = {}

    with open(filepath, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            try:
                cust_id = int(row["客戶"])
                x = float(row["x"])
                y = float(row["y"])
                coords[cust_id] = (x, y)
            except:
                continue

    # 加 depot
    coords[0] = (50.0, 25.0)

    return coords

model = Model("ChainOsr")
Z = 4    # 區域數量
N = 12  # 客戶數量
N0 = 49 # 節點數量(包含倉庫和客戶節點)
SCALE = 1
O = 6   # 旋轉方向
L, W, H = 300,170,165  #貨櫃長寬高
M = 1e6  # 極大值
C0 = 6   #每單位才積委外的費用
C1 = 5.41
P0 = 1000
vi = {}  #客戶 i 所有貨物的總才積
vit = {} #客戶 i 的第 t 件貨物的才積

lito = {} #客戶 i 的第 t 件貨物以擺放方向 o 時的長
wito = {} #客戶 i 的第 t 件貨物以擺放方向 o 時的寬
hito = {} #客戶 i 的第 t 件貨物以擺放方向 o 時的高
pito = {} #擺放方向指標
fit = {}  #脆弱性貨物指標
Gi = {}   #客戶 i 的貨物數量
coord_path = r"datasets/N12_A4_S20250102/generatedCoordinates.csv"
coords = read_generated_coords(coord_path)

print(coords[0])   # depot
print(coords[1])   # customer 1
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
        totalVolume = float(row['總才積'].strip())
        Gi[customer_id] = cargo_count
        vi[customer_id] = totalVolume * SCALE

with open(os.path.join(base_dir, "goods.csv"), newline='', encoding='utf-8-sig') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        customer_id = int(row['客戶'])
        t = int(row['貨物'])
        cargoVolume = float(row['材積'])

        good_l = float(row['長'])
        good_w = float(row['寬'])
        good_h = float(row['高'])
        fragile = int(row['脆弱性'])  # 0 或 1

        volume_scale = SCALE   # 例如 1.0 / 1.25 / 1.5 / 1.75 / 2.0
        length_scale = volume_scale ** (1/3)

        good_l *= length_scale
        good_w *= length_scale
        good_h *= length_scale
        cargoVolume *= volume_scale

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
with open(os.path.join(base_dir, "routes_star.csv"), newline='', encoding='utf-8-sig') as csvfile:
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
with open(os.path.join(base_dir, "routeArcs_star.csv"), newline='', encoding='utf-8-sig') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        area = int(row['區域'])
        start = int(row['起點'])
        end = int(row['終點'])
        Ab[area+1].append((start, end))  
    Ab = OrderedDict(sorted(Ab.items())) 
# print(Ab)

# g as parameter, indexed by (b,i,j)

g = {}
for b in Ab:
    for (i, j) in Ab[b]:
        xi, yi = coords[i]
        xj, yj = coords[j]
        g[b, i, j] = math.hypot(xi - xj, yi - yj)

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

si = model.addVars(
    range(1, N + 1),
    vtype=GRB.BINARY,
    name="si"
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


# --- New variable: epsilon_b ---
eps = model.addVars(range(1, Z+1), vtype=GRB.BINARY, name="eps")
# eps[b] = 1 if area b has any self-owned served customer

# --- New variable: h_b (whether self-owned truck in area b is used) ---
hb = model.addVars(range(1, Z+1), vtype=GRB.BINARY, name="hb")

# overflow indicators
ux, uy, uz = {}, {}, {}
for i in range(1, N + 1):
    for t in range(1, Gi[i] + 1):
        ux[i, t] = model.addVar(vtype=GRB.BINARY, name=f"ux[{i},{t}]")
        uy[i, t] = model.addVar(vtype=GRB.BINARY, name=f"uy[{i},{t}]")
        uz[i, t] = model.addVar(vtype=GRB.BINARY, name=f"uz[{i},{t}]")




outsourceFee = {}
for i in range(1, N + 1):
    vol = vi[i]              # 客戶 i 的總才積
    v = math.ceil(vol)       # 不到 1 才以 1 才計
    outsourceFee[i] = 100 + 30 * max(0, v - 1)
#目標式
obj1 = quicksum(outsourceFee[i] * omega[i] for i in range(1, N+1))   # 第一層目標（較重要）
obj1 += quicksum(
    g[b, i, j] * psi[i, j, b] *C1 
    for b in Ab              # 用 Ab 的 key
    for (i, j) in Ab[b]
)

obj1 += quicksum(P0 * omega[i] for i in range(1, N+1)) 
  

model.setObjective(obj1, GRB.MINIMIZE)

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
        model.addConstr(delta[i, b] == 1)
#3C
for b in range(2, Z + 1):  # Z2 ~ Zu
    if (b-1) in Nobb1:
        for i in Nobb1[b-1]:
            model.addConstr(delta[i, b] == zeta[i, b])
            model.addConstr(delta[i, b-1] == zeta[i, b-1])
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
                psi[i, j, b] <= (delta[i, b] - kappa[i, b]  + delta[j, b] - kappa[j, b]) / 2,
            )
#3G
for b in Ab:
    # 建立 incoming mapping：j → 所有 (i, j)
    incoming = defaultdict(list)
    for (i, j) in Ab[b]:
        incoming[j].append((i, j))

    for j in Nb[b]:  # j ∈ N_b
        model.addConstr(
            quicksum(psi[i, j, b] for (i, j) in incoming.get(j, [])) == delta[j, b] - kappa[j, b],
        )

#3H
for b in Ab:
    # 建立 incoming mapping：j → 所有 (i, j)
    incoming = defaultdict(list)
    for (j, k) in Ab[b]:
        incoming[j].append((j, k))

    for j in Nb[b]:  # j ∈ N_b
        model.addConstr(
            quicksum(psi[j, k, b] for (j, k) in incoming.get(j, [])) == delta[j, b] - kappa[j, b],
        )
# (3I)  delta[j,b] - kappa[j,b] <= eps[b]
for b in Ab:
    for j in Nb.get(b, set()):
        if 1 <= j <= N:  # 只對客戶節點，避開 depot=0 或其他非客戶編號
            model.addConstr(delta[j, b] - kappa[j, b] <= eps[b],
                            name=f"3I[{j},{b}]")

# (3J)  sum_{j in Nb[b]} (delta[j,b] - kappa[j,b]) >= eps[b]
for b in Ab:
    model.addConstr(
        quicksum(delta[j, b] - kappa[j, b]
                    for j in Nb.get(b, set())
                    if 1 <= j <= N) >= eps[b],
        name=f"3J[{b}]"
    )
# (3K1) depot outgoing flow equals h[b]
for b in Ab:
    model.addConstr(
        quicksum(psi[0, v, b] for (u, v) in Ab.get(b, set()) if u == 0) == hb[b],
        name=f"3K1[{b}]"
    )


# (3K2) depot incoming flow equals h[b]
for b in Ab:
    model.addConstr(
        quicksum(psi[u, 0, b] for (u, v) in Ab.get(b, set()) if v == 0) == hb[b],
        name=f"3K2[{b}]"
    )

for b in Ab:
    model.addConstr(hb[b] == eps[b], name=f"link_hb_eps[{b}]")

#4
for (i, j) in A:
    if i != 0 and j != 0:
        model.addConstr(
            quicksum(psi[i, j, b] for b in Ab if (i, j) in Ab[b]) <= gamma[i, j]
        )

#gamma設計
for i in range(1, N+1): 
    for j in range(1, N+1): 
        for k in range(1, N+1): 
            if i != j and i != k and j != k: 
                model.addConstr(
                    gamma[i, j] + gamma[j, k] <= 1 + gamma[i, k],
                )

for i in range(1, N+1): 
    for j in range(1, N+1):
        if i != j:
           model.addConstr(gamma[i, j] + gamma[j, i] <= 1) 

for i in range(1, N+1):  
    for j in range(1, N+1):  
        if i != j:  
            for b in range(1, Z + 1):  
                for bprime in range(1, Z + 1):
                    if b != bprime: 
                        model.addConstr(
                            gamma[i, j] + gamma[j, i] <= 2 - delta[i, b] - delta[j, bprime],
                        )

for i in range(1, N+1):  
    for j in range(1, N+1):  
        if i != j:  
            for b in range(1, Z + 1):  
                model.addConstr(
                    gamma[i, j] + gamma[j, i] >= 1 - (2 - delta[i, b] - delta[j, b])
                )
#5A
for i in range(1, N+1):  
    for j in range(1, N+1):  
        if i != j:
            for t in range(1, Gi[i] + 1):
                for tprime in range(1, Gi[j] + 1):
                    model.addConstr(
                        xprime[i, t, j, tprime] + xprime[j, tprime, i, t] +
                        yprime[i, t, j, tprime] + yprime[j, tprime, i, t] + 
                        zprime[i, t, j, tprime] + zprime[j, tprime, i, t] <= 3 * (gamma[i, j] + gamma[j, i]),
                    )

#5B
for i in range(1, N+1):  
    for j in range(1, N+1):  
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

#6B~D (修正版：外包時退出 non-overlap)
for i in range(1, N+1):
    for t in range(1, Gi[i] + 1):
        for j in range(1, N+1):
            for tprime in range(1, Gi[j] + 1):
                if (i == j and t != tprime) or (i != j):

                    size_i_L = quicksum(alpha[i, t, o] * lito[i][t-1][o] for o in range(O))
                    size_i_W = quicksum(alpha[i, t, o] * wito[i][t-1][o] for o in range(O))
                    size_i_H = quicksum(alpha[i, t, o] * hito[i][t-1][o] for o in range(O))

                    model.addConstr(
                        x[i, t] + size_i_L
                        <= x[j, tprime] + M * (1 - xprime[i, t, j, tprime]) 
                    )
                    model.addConstr(
                        y[i, t] + size_i_W
                        <= y[j, tprime] + M * (1 - yprime[i, t, j, tprime]) 
                    )
                    model.addConstr(
                        z[i, t] + size_i_H
                        <= z[j, tprime] + M * (1 - zprime[i, t, j, tprime]) 
                    )

    
#6E~F 
for i in range(1, N + 1):
    for t in range(1, Gi[i] + 1):
        model.addConstr(
            x[i, t] + quicksum(alpha[i, t, o] * lito[i][t-1][o] for o in range(O)) <= L + M*omega[i]
        )
        model.addConstr(
            y[i, t] + quicksum(alpha[i, t, o] * wito[i][t-1][o] for o in range(O)) <= W  + M*omega[i]
        )
        model.addConstr(
            z[i, t] + quicksum(alpha[i, t, o] * hito[i][t-1][o] for o in range(O)) <= H  + M*omega[i]
        )


#26P
#for (i, j) in A:
 #   if i != 0 and j != 0:
  #      model.addConstr(omega[j] >= omega[i] + gamma[i,j] - 1)
for b in range(1, Z+1):
    for i in range(1, N+1):
        for j in range(1, N+1):
            if i != j:
                model.addConstr(
                    kappa[j,b] >= kappa[i,b] + gamma[i,j] + delta[i,b] + delta[j,b] - 3
                )



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
                
print("\n==============================")
print(" Solution Check Prints (26/3/4/5)")
print("==============================")

# 0) 狀態與目標值
status = model.Status
print("Model status =", status)
if status == GRB.OPTIMAL:
    print("Objective =", model.ObjVal)
elif status in (GRB.INFEASIBLE, GRB.INF_OR_UNBD):
    print("Infeasible. (Suggestion) compute IIS:")
    print("  model.computeIIS(); model.write('case.ilp')")
    # 你也可以直接打開這兩行跑 IIS：
    # model.computeIIS()
    # model.write("case.ilp")
    raise SystemExit

# ------------------------------------------------
# (26) omega：是否外包/超長放鬆（你要看哪些客戶被拉成 1）
# ------------------------------------------------
print("\n(26) omega[i] (outsourcing / relax length):")
omega_ones = []
try:
    for i in omega.keys():
        if omega[i].X > 0.5:
            omega_ones.append(i)
    print("omega=1 customers:", omega_ones if omega_ones else "None")
except Exception as e:
    print("omega not found / name mismatch:", e)

# ------------------------------------------------
# (3) delta：客戶指派到哪個區域/車
#    你要看：每個客戶是否剛好指派一次、以及是否落在允許服務區域
# ------------------------------------------------
print("\n(3) delta[i,b] assignment:")
try:
    # 推出 customers 清單：用 delta 的 key 或 Gi 的 key 都行
    # 這裡用 delta 的 key 最直觀
    customers = sorted({i for (i, b) in delta.keys()})
    areas = sorted({b for (i, b) in delta.keys()})

    # 每個客戶被指派到哪個區域
    assigned = {}
    bad_multi = []
    bad_none = []
    for i in customers:
        chosen = [b for b in areas if (i, b) in delta and delta[i, b].X > 0.5]
        if len(chosen) == 1:
            assigned[i] = chosen[0]
        elif len(chosen) == 0:
            bad_none.append(i)
        else:
            bad_multi.append((i, chosen))

    print("Assigned area per customer (i -> b):")
    for i in customers:
        print(f"  {i} -> {assigned.get(i, None)}")

    if bad_none:
        print("  [FAIL] customers with no assigned area:", bad_none)
    if bad_multi:
        print("  [FAIL] customers assigned to multiple areas:", bad_multi)
    if (not bad_none) and (not bad_multi):
        print("  PASS: each customer assigned to exactly one area.")

    # 允許服務區域檢核（如果你有 serviceAreas / allowedAreas）
    # 你程式裡若有類似 allowedAreas[i] = set([...])，就可以直接檢查。
    # 這裡提供通用寫法：如果你有 service_area_flags[i] = [0/1,...] 也能改。
    if "allowedAreas" in globals():
        bad_allowed = []
        for i, b in assigned.items():
            if b not in allowedAreas[i]:
                bad_allowed.append((i, b, sorted(list(allowedAreas[i]))))
        if bad_allowed:
            print("  [FAIL] assigned area not allowed:")
            for i, b, allow in bad_allowed:
                print(f"    customer {i} assigned {b}, allowed={allow}")
        else:
            print("  PASS: assigned area within allowed set.")
    else:
        print("  NOTE: allowedAreas not found in code; skip 'allowed service area' check.")
except Exception as e:
    print("delta assignment check failed (name mismatch / key format):", e)

# ------------------------------------------------
# (4) psi & gamma：路由弧與先後關係
#    你要看：選了哪些弧、gamma 是否變成一致的順序（同區域二選一）
# ------------------------------------------------
print("\n(4) psi[i,j,b] selected arcs (only print X=1):")
try:
    cnt = 0
    for key in psi.keys():
        if psi[key].X > 0.5:
            print("  psi", key, "=1")
            cnt += 1
    if cnt == 0:
        print("  (none)")
except Exception as e:
    print("psi not found / name mismatch:", e)

print("\n(4/5) gamma[i,j] precedence (only print X=1 for i<j):")
try:
    cnt = 0
    # 只印 i<j 避免重複
    for (i, j) in gamma.keys():
        if i < j and gamma[i, j].X > 0.5:
            print(f"  gamma[{i},{j}] = 1")
            cnt += 1
    if cnt == 0:
        print("  (none or all 0)")
except Exception as e:
    print("gamma not found / name mismatch:", e)

# ------------------------------------------------
# (5) 裝載：座標 x,y,z + 邊界檢查（粗略）
#    你要看：是否有跑出位置、是否超出車廂、以及懸浮(模型通常允許)
# ------------------------------------------------
print("\n(5) packing coordinates & bounds check:")
L, W, H = 300, 170, 165  # 你車廂尺寸
out_of_box = []
try:
    # 你 x/y/z 的 key 可能是 (i,t) 或 (i,t,v) 之類
    # 這裡假設至少有 (i,t) 這種維度
    for key in x.keys():
        xv = x[key].X
        yv = y[key].X
        zv = z[key].X
        if xv < -1e-6 or yv < -1e-6 or zv < -1e-6 or xv > L + 1e-6 or yv > W + 1e-6 or zv > H + 1e-6:
            out_of_box.append((key, xv, yv, zv))
    if out_of_box:
        print("  [FAIL] some items out of container bounds:")
        for key, xv, yv, zv in out_of_box[:50]:
            print(f"    item {key} -> ({xv:.2f}, {yv:.2f}, {zv:.2f})")
    else:
        print("  PASS: all (x,y,z) are within container range (rough).")
except Exception as e:
    print("x/y/z not found / key mismatch:", e)

print("\nNOTE:")
print(" - 目前只做『座標落點』的粗略邊界檢查；要完整驗證不重疊與 x+placed_L<=L 等，必須輸出每件貨物的旋轉後尺寸(placed_L/W/H)。")
print("==============================\n")

# print("\n=== 各區貨物材積總和")
# for b in range(1, Z + 1):
#     if R[b].X > 0.01:  # 避免浮點誤差影響
#         print(f"區域 {b} 的總裝載材積為：{R[b].X:.2f} 單位³")
#     else:
#         print(f"區域 {b} 的總裝載材積為：0 單位³")

print("\n==============================")
print(" Obj1 Breakdown Check ")
print("==============================")

# --- (A1) Outsourcing Base Cost ---
outsourcing_base_cost = 0
outsourcing_base_detail = []

for i in range(1, N + 1):
    if omega[i].X > 0.5:
        cost_i = outsourceFee[i]
        outsourcing_base_cost += cost_i
        outsourcing_base_detail.append((i, vi[i], outsourceFee[i], cost_i))

print("\n[Obj1-A1] Outsourcing Base Cost =", outsourcing_base_cost)
print("[Obj1-A1 Details] (customer, vi, outsourceFee[i], contribution)")
for item in outsourcing_base_detail:
    print(" ", item)


# --- (A2) Outsourcing Penalty Cost ---
penalty_cost = 0
penalty_detail = []

for i in range(1, N + 1):
    if omega[i].X > 0.5:
        cost_i = P0
        penalty_cost += cost_i
        penalty_detail.append((i, cost_i))

print("\n[Obj1-A2] Outsourcing Penalty Cost =", penalty_cost)
print("[Obj1-A2 Details] (customer, penalty)")
for item in penalty_detail:
    print(" ", item)


# --- (B) Self-owned Fuel / Arc Cost ---
fuel_cost = 0
selected_arcs = []

for b in Ab:
    for (i, j) in Ab[b]:
        if psi[i, j, b].X > 0.5:
            arc_cost = C1 * g[b, i, j]
            fuel_cost += arc_cost
            selected_arcs.append((b, i, j, g[b, i, j], arc_cost))

print("\n[Obj1-B] Fuel / Arc Cost =", fuel_cost)
print("[Obj1-B Details] (b, i, j, distance, contribution)")
for arc in selected_arcs:
    print(" ", arc)


# --- Final Check ---
total_manual = outsourcing_base_cost + penalty_cost + fuel_cost

print("\n[Obj1 Manual Sum] =", total_manual)

try:
    print("[Obj1 Model Value] =", obj1.getValue())
except:
    print("[Obj1 Model Value] cannot access obj1 directly (multi-objective mode)")


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

# ---- Debug: print delta/kappa/served for outsourced customers ----
def v(x):
    """Safe value getter for Gurobi vars (handles None)."""
    if x is None:
        return None
    try:
        return x.X
    except:
        try:
            return x.x
        except:
            return None

# 1) 找出外包客戶
outsourced = []
for i in range(1, N + 1):
    if v(omega[i]) is not None and v(omega[i]) > 0.5:
        outsourced.append(i)

print("\n[DEBUG] Outsourced customers (omega=1):", outsourced if outsourced else "None")

# 2) 對每個外包客戶，印 delta/kappa/delta-kappa
#    b 的集合用你模型裡實際存在的區域 key：優先用 Ab.keys()，不然用 1..Z
try:
    B_list = sorted(list(Ab.keys()))
except:
    B_list = list(range(1, Z + 1))

for i in outsourced:
    print(f"\n[DEBUG] Customer i={i}, omega={v(omega[i])}")

    # 找 i 被分到的區（delta=1 的那個）
    assigned_bs = []
    for b in B_list:
        if (i, b) in delta:
            dv = v(delta[i, b])
            if dv is not None and dv > 0.5:
                assigned_bs.append(b)
    print("  assigned b (delta=1):", assigned_bs if assigned_bs else "None")

    # 印每個 b 的 delta/kappa/served
    for b in B_list:
        if (i, b) not in delta or (i, b) not in kappa:
            continue
        dv = v(delta[i, b])
        kv = v(kappa[i, b])
        if dv is None or kv is None:
            continue
        served = dv - kv
        # 只印有關聯的（delta 或 kappa 有 1 的）
        if dv > 0.5 or kv > 0.5:
            print(f"  b={b}: delta={dv:.0f}, kappa={kv:.0f}, delta-kappa={served:.0f}")

area_stacks = {b: [] for b in range(1, Z + 1)}
    
i = 12
print(f"\n[DEBUG customer {i}]")
print("omega =", omega[i].X)
print("si    =", si[i].X)

for t in range(1, Gi[i] + 1):
    placed_L = sum(alpha[i, t, o].X * lito[i][t-1][o] for o in range(O))
    placed_W = sum(alpha[i, t, o].X * wito[i][t-1][o] for o in range(O))
    placed_H = sum(alpha[i, t, o].X * hito[i][t-1][o] for o in range(O))

    print(f"  item t={t}")
    print("    x,y,z =", x[i, t].X, y[i, t].X, z[i, t].X)
    print("    placed_L/W/H =", placed_L, placed_W, placed_H)

    if (i, t) in ux:
        print("    ux,uy,uz =", ux[i, t].X, uy[i, t].X, uz[i, t].X)

    print("    x+L - L =", x[i, t].X + placed_L - L)
    print("    y+W - W =", y[i, t].X + placed_W - W)
    print("    z+H - H =", z[i, t].X + placed_H - H)

    print("\n[DEBUG outsourced customers]")
for i in range(1, N + 1):
    if omega[i].X > 0.5:
        print(f"\ncustomer {i}: omega={omega[i].X}, si={si[i].X}")
        for t in range(1, Gi[i] + 1):
            vals = []
            if (i, t) in ux:
                vals = [ux[i, t].X, uy[i, t].X, uz[i, t].X]
            print(f"  t={t}, overflow flags={vals}")

    
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