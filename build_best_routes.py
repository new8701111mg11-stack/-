from gurobipy import Model, GRB, quicksum
from collections import defaultdict, OrderedDict
import csv
import math
import os
from typing import Dict, List, Tuple, Set

# ------------------------------------------------------------
# 讀 Solomon 座標檔
# ------------------------------------------------------------
def read_solomon_xy(filepath: str) -> Dict[int, Tuple[float, float]]:
    coords = {}
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    start = None
    for i, line in enumerate(lines):
        if "CUST NO." in line and "XCOORD" in line:
            start = i + 1
            break
    if start is None:
        raise ValueError("找不到 CUSTOMER 表頭")

    for line in lines[start:]:
        parts = line.strip().split()
        if len(parts) < 3:
            continue
        try:
            cust_id = int(parts[0])
            x = float(parts[1])
            y = float(parts[2])
        except Exception:
            continue
        coords[cust_id] = (x, y)
    return coords


# ------------------------------------------------------------
# 讀你的既有資料：customerInfo / serviceArea
# 這支程式保留你原本的資料命名風格，方便接回主模型
# ------------------------------------------------------------
def read_instance(base_dir: str, solomon_path: str, Z: int):
    coords = read_solomon_xy(solomon_path)

    vi = {}
    Gi = {}
    Neb = defaultdict(set)      # 非重疊區域客戶
    Nobb1 = defaultdict(list)   # b 與 b+1 重疊客戶

    with open(os.path.join(base_dir, "customerInfo.csv"), newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            customer_id = int(row['客戶'].strip())
            cargo_count = int(row['貨物數'].strip())
            total_volume = float(row['總才積'].strip())
            Gi[customer_id] = cargo_count
            vi[customer_id] = total_volume

    with open(os.path.join(base_dir, "serviceArea.csv"), newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            customer_id = int(row['客戶'])
            service_vector = [int(row[f'服務區域{b}'].strip()) for b in range(1, Z + 1)]

            if service_vector.count(1) == 1:
                area_index = service_vector.index(1) + 1
                Neb[area_index].add(customer_id)

            for b in range(Z - 1):
                if service_vector[b] == 1 and service_vector[b + 1] == 1:
                    Nobb1[b + 1].append(customer_id)

    Neb = OrderedDict((k, sorted(v)) for k, v in sorted(Neb.items()))
    Nobb1 = OrderedDict((k, sorted(v)) for k, v in sorted(Nobb1.items()))

    return coords, vi, Gi, Neb, Nobb1


# ------------------------------------------------------------
# 對每個區域 r 建立「包含所有可能 OSR 客戶」的集合
# 依你老師文件的意思：
#   r 的客戶 = Neb[r] + 左側 overlap + 右側 overlap
# ------------------------------------------------------------
def build_vrp_customer_sets(Z: int,
                            Neb: Dict[int, List[int]],
                            Nobb1: Dict[int, List[int]]) -> Dict[int, List[int]]:
    area_customers = {}
    for r in range(1, Z + 1):
        s = set(Neb.get(r, []))

        # 左側重疊：r-1 與 r 的 overlap 存在 Nobb1[r-1]
        if r - 1 in Nobb1:
            s.update(Nobb1[r - 1])

        # 右側重疊：r 與 r+1 的 overlap 存在 Nobb1[r]
        if r in Nobb1:
            s.update(Nobb1[r])

        area_customers[r] = sorted(s)
    return area_customers


# ------------------------------------------------------------
# 距離
# ------------------------------------------------------------
def euclidean(coords: Dict[int, Tuple[float, float]], i: int, j: int) -> float:
    xi, yi = coords[i]
    xj, yj = coords[j]
    return math.hypot(xi - xj, yi - yj)


# ------------------------------------------------------------
# 單一區域的 TSP/單車 VRP
# - 1 台車
# - depot = 0
# - 客戶都要拜訪一次
# - 回到 depot
# - 不考慮容量、不考慮 3D 裝載
# 這正符合你老師文件裡「先產生 R*_r」的用途
# ------------------------------------------------------------
def solve_single_area_tsp(area_id: int,
                          customers: List[int],
                          coords: Dict[int, Tuple[float, float]],
                          time_limit: int = 60,
                          mip_gap: float = 0.0):
    model = Model(f"Area_{area_id}_BestRoute")
    model.Params.OutputFlag = 0
    model.Params.TimeLimit = time_limit
    model.Params.MIPGap = mip_gap

    if not customers:
        return {
            "area": area_id,
            "customers": [],
            "objective": 0.0,
            "arcs": [],
            "route": [0, 0],
            "status": "EMPTY"
        }

    nodes = [0] + customers
    arcs = [(i, j) for i in nodes for j in nodes if i != j]
    d = {(i, j): euclidean(coords, i, j) for (i, j) in arcs}

    x = model.addVars(arcs, vtype=GRB.BINARY, name="x")

    # MTZ 順序變數，只對 customer 定義
    u = model.addVars(customers, lb=1, ub=len(customers), vtype=GRB.CONTINUOUS, name="u")

    model.setObjective(quicksum(d[i, j] * x[i, j] for (i, j) in arcs), GRB.MINIMIZE)

    # 每個客戶 exactly one in / one out
    for j in customers:
        model.addConstr(quicksum(x[i, j] for i in nodes if i != j) == 1, name=f"in_{j}")
        model.addConstr(quicksum(x[j, k] for k in nodes if k != j) == 1, name=f"out_{j}")

    # depot exactly one departure / one return
    model.addConstr(quicksum(x[0, j] for j in customers) == 1, name="depot_out")
    model.addConstr(quicksum(x[i, 0] for i in customers) == 1, name="depot_in")

    # MTZ subtour elimination
    n = len(customers)
    if n >= 2:
        for i in customers:
            for j in customers:
                if i != j:
                    model.addConstr(u[i] - u[j] + n * x[i, j] <= n - 1,
                                    name=f"mtz_{i}_{j}")

    model.optimize()

    if model.Status not in [GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.SUBOPTIMAL]:
        return {
            "area": area_id,
            "customers": customers,
            "objective": None,
            "arcs": [],
            "route": None,
            "status": f"STATUS_{model.Status}"
        }

    selected_arcs = [(i, j) for (i, j) in arcs if x[i, j].X > 0.5]
    route = extract_route_from_arcs(selected_arcs)

    return {
        "area": area_id,
        "customers": customers,
        "objective": model.ObjVal,
        "arcs": selected_arcs,
        "route": route,
        "status": "OK"
    }


# ------------------------------------------------------------
# 從弧集合重建路徑
# 輸出例: [0, 7, 6, 12, 16, 14, 11, 8, 13, 15, 0]
# ------------------------------------------------------------
def extract_route_from_arcs(arcs: List[Tuple[int, int]]) -> List[int]:
    if not arcs:
        return []

    succ = {i: j for i, j in arcs}
    route = [0]
    visited = set([0])
    cur = 0

    while True:
        nxt = succ.get(cur)
        if nxt is None:
            break
        route.append(nxt)
        if nxt == 0:
            break
        if nxt in visited:
            # 保底，避免異常循環
            break
        visited.add(nxt)
        cur = nxt

    return route


# ------------------------------------------------------------
# 轉成你比較好看的 R*_r 格式
# 例: (7-6)-(12)-(16-14)-(11)-(8)-(13)-(15)
# 目前先簡單輸出成 0-a-b-c-0 的字串與客戶序列
# 若你之後要完全模仿老師圖上的分組，我再幫你改
# ------------------------------------------------------------
def route_to_string(route: List[int]) -> str:
    if not route:
        return ""
    return "-".join(map(str, route))


# ------------------------------------------------------------
# 主程式
# ------------------------------------------------------------
def generate_forward_arcs_from_route(route: List[int]) -> List[Tuple[int, int]]:
    """
    例如 route = [0, 6, 7, 12, 0]
    產生:
    (0,6), (0,7), (0,12),
    (6,7), (6,12), (6,0),
    (7,12), (7,0),
    (12,0)
    """
    arcs = []
    n = len(route)
    for i in range(n - 1):
        for j in range(i + 1, n):
            arcs.append((route[i], route[j]))
    return arcs


def write_routes_csv(routes_by_zone: Dict[int, List[int]], out_path: str):
    """
    輸出成和原本 routes.csv 一樣的格式
    """
    max_len = max(len(route) for route in routes_by_zone.values())

    header = ["區域"] + [f"節點{i}" for i in range(1, max_len + 1)]

    with open(out_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(header)

        for zone in sorted(routes_by_zone.keys()):
            route = routes_by_zone[zone]

            # 原本檔案區域是 0-based，所以輸出 zone-1
            row = [zone - 1] + route + [""] * (max_len - len(route))
            writer.writerow(row)


def write_route_arcs_csv(routes_by_zone: Dict[int, List[int]], out_path: str):
    """
    輸出成和原本 routeArcs.csv 一樣的格式
    """
    with open(out_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["區域", "起點", "終點"])

        for zone in sorted(routes_by_zone.keys()):
            route = routes_by_zone[zone]
            arcs = generate_forward_arcs_from_route(route)

            for i, j in arcs:
                # 原本檔案區域是 0-based，所以輸出 zone-1
                writer.writerow([zone - 1, i, j])


def main():
    # ===== 你可以改這裡 =====
    Z = 4
    base_dir = r"datasets/N12_A4_S20250102"
    solomon_path = r"datasets/C101.txt"

    output_txt = os.path.join(base_dir, "best_routes_star.txt")
    output_routes_csv = os.path.join(base_dir, "routes_star.csv")
    output_route_arcs_csv = os.path.join(base_dir, "routeArcs_star.csv")
    # ======================

    coords, vi, Gi, Neb, Nobb1 = read_instance(base_dir, solomon_path, Z)
    area_customers = build_vrp_customer_sets(Z, Neb, Nobb1)

    print("=== 各區『包含所有可能 OSR 客戶』集合 ===")
    for r in range(1, Z + 1):
        print(f"R{r}: {area_customers[r]}")

    results = []
    for r in range(1, Z + 1):
        res = solve_single_area_tsp(
            area_id=r,
            customers=area_customers[r],
            coords=coords,
            time_limit=60,
            mip_gap=0.0,
        )
        results.append(res)

    print("\n=== 各區最佳路徑 R*_r ===")
    for res in results:
        print(f"區域 {res['area']}")
        print(f"  customers : {res['customers']}")
        print(f"  route     : {res['route']}")
        print(f"  route str : {route_to_string(res['route'])}")
        print(f"  obj dist  : {res['objective']}")
        print(f"  status    : {res['status']}")

    # 輸出 txt
    with open(output_txt, "w", encoding="utf-8") as f:
        f.write("=== 各區最佳路徑 R*_r ===\n")
        for res in results:
            f.write(f"區域 {res['area']}\n")
            f.write(f"customers : {res['customers']}\n")
            f.write(f"route     : {res['route']}\n")
            f.write(f"route str : {route_to_string(res['route'])}\n")
            f.write(f"obj dist  : {res['objective']}\n")
            f.write(f"status    : {res['status']}\n\n")

    # 整理成 {區域: 路徑}，給 routes_star.csv / routeArcs_star.csv 用
    best_routes = {}
    for res in results:
        best_routes[res["area"]] = res["route"]

    # 輸出成和原本一樣格式的 CSV
    write_routes_csv(best_routes, output_routes_csv)
    write_route_arcs_csv(best_routes, output_route_arcs_csv)

    print(f"\n已輸出：{output_txt}")
    print(f"已輸出：{output_routes_csv}")
    print(f"已輸出：{output_route_arcs_csv}")


if __name__ == "__main__":
    main()