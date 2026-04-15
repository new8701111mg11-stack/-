import os
import csv
import math
from gurobipy import Model, GRB, quicksum


# =========================================================
# 基本路徑設定
# =========================================================
base_dir = r"datasets/N12_A4_S20250102"

case_file = "saved_sol_3_case_0_area3_truck3_tryCust3.txt"
coord_file = "generatedCoordinates.csv"

case_path = os.path.join(base_dir, case_file)
coord_path = os.path.join(base_dir, coord_file)


# =========================================================
# 1. 讀 generatedCoordinates.csv
# =========================================================
def read_generated_coords(base_dir, filename="generatedCoordinates.csv"):
    filepath = os.path.join(base_dir, filename)
    coords = {}

    with open(filepath, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            try:
                cust_id = int(row["客戶"].strip())
                x = float(row["x"].strip())
                y = float(row["y"].strip())
                coords[cust_id] = (x, y)
            except:
                continue

    # depot
    coords[0] = (50.0, 25.0)
    return coords


# =========================================================
# 2. 讀單一 case txt
# =========================================================
def read_single_truck_case(base_dir, filename):
    filepath = os.path.join(base_dir, filename)

    data = {
        "case_id": None,
        "base_fitness": None,
        "area": None,
        "truck_id": None,
        "trying_customer": None,
        "failed_customer": None,
        "failed_cargo": None,
        "container": None,
        "candidate_count": 0,
        "boundary_fail_count": 0,
        "collision_fail_count": 0,
        "best_boundary_violation": 0.0,
        "best_total_overlap": 0.0,
        "collision_dominant": 0,

        "loaded_customers_before_fail": [],
        "truck_case_customers": [],
        "full_planned_customers": [],
        "remaining_customers_after_fail": [],
        "fail_index_in_planned_route": -1,

        "items": []
    }

    in_items = False

    with open(filepath, "r", encoding="utf-8") as fin:
        for raw in fin:
            line = raw.strip()
            if not line:
                continue

            if line == "ITEMS_BEGIN":
                in_items = True
                continue

            if line == "ITEMS_END":
                in_items = False
                continue

            if in_items:
                if line.startswith("#"):
                    continue

                parts = line.split()

                # 新格式：
                # customerId cargoId l w h isPreloaded itemRole ori1 ori2 ori3 ori4 ori5 ori6
                if len(parts) >= 13:
                    item = {
                        "customer_id": int(parts[0]),
                        "cargo_id": int(parts[1]),
                        "l": int(parts[2]),
                        "w": int(parts[3]),
                        "h": int(parts[4]),
                        "is_preloaded": int(parts[5]),
                        "item_role": int(parts[6]),
                        "orientation_flags": [
                            int(parts[7]), int(parts[8]), int(parts[9]),
                            int(parts[10]), int(parts[11]), int(parts[12])
                        ]
                    }
                    data["items"].append(item)
                continue

            parts = line.split()
            key = parts[0]

            if key == "CASE_ID":
                data["case_id"] = int(parts[1])

            elif key == "FITNESS":
                data["base_fitness"] = float(parts[1])

            elif key == "BASE_FITNESS":
                data["base_fitness"] = float(parts[1])

            elif key == "AREA":
                data["area"] = int(parts[1])

            elif key == "TRUCK_ID":
                data["truck_id"] = int(parts[1])

            elif key == "TRYING_CUSTOMER":
                data["trying_customer"] = int(parts[1])

            elif key == "FAILED_CUSTOMER":
                data["failed_customer"] = int(parts[1])

            elif key == "FAILED_CARGO":
                data["failed_cargo"] = int(parts[1])

            elif key == "CONTAINER":
                data["container"] = (int(parts[1]), int(parts[2]), int(parts[3]))

            elif key == "CANDIDATE_COUNT":
                data["candidate_count"] = int(parts[1])

            elif key == "BOUNDARY_FAIL_COUNT":
                data["boundary_fail_count"] = int(parts[1])

            elif key == "COLLISION_FAIL_COUNT":
                data["collision_fail_count"] = int(parts[1])

            elif key == "BEST_BOUNDARY_VIOLATION":
                data["best_boundary_violation"] = float(parts[1])

            elif key == "BEST_TOTAL_OVERLAP":
                data["best_total_overlap"] = float(parts[1])

            elif key == "COLLISION_DOMINANT":
                data["collision_dominant"] = int(parts[1])

            elif key == "LOADED_CUSTOMERS_BEFORE_FAIL":
                data["loaded_customers_before_fail"] = list(map(int, parts[1:]))

            elif key == "TRUCK_CASE_CUSTOMERS":
                data["truck_case_customers"] = list(map(int, parts[1:]))

            elif key == "FULL_PLANNED_CUSTOMERS":
                data["full_planned_customers"] = list(map(int, parts[1:]))

            elif key == "REMAINING_CUSTOMERS_AFTER_FAIL":
                data["remaining_customers_after_fail"] = list(map(int, parts[1:]))

            elif key == "FAIL_INDEX_IN_PLANNED_ROUTE":
                data["fail_index_in_planned_route"] = int(parts[1])

    return data

# =========================================================
# 3. orientation 尺寸
# =========================================================
def get_oriented_dims(l, w, h):
    return {
        0: (l, w, h),
        1: (l, h, w),
        2: (w, l, h),
        3: (w, h, l),
        4: (h, l, w),
        5: (h, w, l),
    }


# =========================================================
# 4. 客戶總才積
#    目前先從 case item 用 l*w*h 累加
#    之後若你要和 C++ 完全一致，可改成另外讀正式 vi
# =========================================================
def build_customer_volumes(case_data):
    vi = {}

    for item in case_data["items"]:
        cid = item["customer_id"]
        vol = item["l"] * item["w"] * item["h"]/27000

        if cid not in vi:
            vi[cid] = 0.0
        vi[cid] += vol

    return vi


# =========================================================
# 5. 外包費
# =========================================================
def build_outsource_fee(vi):
    outsource_fee = {}
    for i, vol in vi.items():
        v = math.ceil(vol)
        outsource_fee[i] = 100 + 30 * max(0, v - 3)
    return outsource_fee


# =========================================================
# 6. 距離
# =========================================================
def euclidean_distance(a, b, coords):
    xa, ya = coords[a]
    xb, yb = coords[b]
    return math.hypot(xa - xb, ya - yb)


def compute_route_distance(route, coords):
    if not route:
        return 0.0

    dist = 0.0
    dist += euclidean_distance(0, route[0], coords)
    for i in range(len(route) - 1):
        dist += euclidean_distance(route[i], route[i + 1], coords)
    dist += euclidean_distance(route[-1], 0, coords)
    return dist


# =========================================================
# 7. 簡化 route
# =========================================================
def build_routes_from_case(case_data):
    old_route = case_data["loaded_customers_before_fail"][:]
    local_rescued_route = case_data["truck_case_customers"][:]
    full_planned_route = case_data["full_planned_customers"][:]
    remaining_after_fail = case_data["remaining_customers_after_fail"][:]
    fail_index = case_data["fail_index_in_planned_route"]

    return old_route, local_rescued_route, full_planned_route, remaining_after_fail, fail_index

def build_single_truck_packing_model(case_data, coords, C1=5.41):
    model = Model(f"single_truck_case_{case_data['case_id']}")

    items = case_data["items"]
    n_items = len(items)

    L, W, H = case_data["container"]

    I = range(n_items)
    R = range(6)

    vi = build_customer_volumes(case_data)
    outsource_fee = build_outsource_fee(vi)

    trying_customer = case_data["trying_customer"]
    base_fitness = case_data["base_fitness"]

    old_route, local_rescued_route, full_planned_route, remaining_after_fail, fail_index = build_routes_from_case(case_data)

    old_dist = compute_route_distance(old_route, coords)
    local_rescued_dist = compute_route_distance(local_rescued_route, coords)
    full_planned_dist = compute_route_distance(full_planned_route, coords)

    # 版本 B：用 full planned route 估完整救回後的額外油耗
    extra_fuel_cost = C1 * max(0.0, full_planned_dist - old_dist)
    saved_outsource_cost = outsource_fee.get(trying_customer, 0.0)

    x = model.addVars(I, vtype=GRB.CONTINUOUS, lb=0.0, name="x")
    y = model.addVars(I, vtype=GRB.CONTINUOUS, lb=0.0, name="y")
    z = model.addVars(I, vtype=GRB.CONTINUOUS, lb=0.0, name="z")

    placed = model.addVars(I, vtype=GRB.BINARY, name="placed")
    o = model.addVars(I, R, vtype=GRB.BINARY, name="o")

    lx = model.addVars(I, vtype=GRB.CONTINUOUS, lb=0.0, name="lx")
    wy = model.addVars(I, vtype=GRB.CONTINUOUS, lb=0.0, name="wy")
    hz = model.addVars(I, vtype=GRB.CONTINUOUS, lb=0.0, name="hz")

    rescue_success = model.addVar(vtype=GRB.BINARY, name="rescue_success")

    for i in I:
        model.addConstr(
            quicksum(o[i, r] for r in R) == placed[i],
            name=f"orient_choose_{i}"
        )

        dims = get_oriented_dims(items[i]["l"], items[i]["w"], items[i]["h"])

        for r in R:
            if items[i]["orientation_flags"][r] == 0:
                model.addConstr(o[i, r] == 0, name=f"forbid_ori_{i}_{r}")

        model.addConstr(
            lx[i] == quicksum(dims[r][0] * o[i, r] for r in R),
            name=f"lx_def_{i}"
        )
        model.addConstr(
            wy[i] == quicksum(dims[r][1] * o[i, r] for r in R),
            name=f"wy_def_{i}"
        )
        model.addConstr(
            hz[i] == quicksum(dims[r][2] * o[i, r] for r in R),
            name=f"hz_def_{i}"
        )

    M = 10**6
    for i in I:
        model.addConstr(x[i] + lx[i] <= L + M * (1 - placed[i]), name=f"bound_x_{i}")
        model.addConstr(y[i] + wy[i] <= W + M * (1 - placed[i]), name=f"bound_y_{i}")
        model.addConstr(z[i] + hz[i] <= H + M * (1 - placed[i]), name=f"bound_z_{i}")

    for i in I:
        model.addConstr(placed[i] == 1, name=f"must_place_{i}")

    model.addConstr(rescue_success == 1, name="force_rescue_success")

    pair_list = [(i, k) for i in I for k in I if i < k]

    left  = model.addVars(pair_list, vtype=GRB.BINARY, name="left")
    right = model.addVars(pair_list, vtype=GRB.BINARY, name="right")
    front = model.addVars(pair_list, vtype=GRB.BINARY, name="front")
    back  = model.addVars(pair_list, vtype=GRB.BINARY, name="back")
    below = model.addVars(pair_list, vtype=GRB.BINARY, name="below")
    above = model.addVars(pair_list, vtype=GRB.BINARY, name="above")

    for (i, k) in pair_list:
        model.addConstr(
            left[i, k] + right[i, k] +
            front[i, k] + back[i, k] +
            below[i, k] + above[i, k]
            >= placed[i] + placed[k] - 1,
            name=f"sep_choose_{i}_{k}"
        )

        model.addConstr(
            x[i] + lx[i] <= x[k] + M * (1 - left[i, k]),
            name=f"left_sep_{i}_{k}"
        )
        model.addConstr(
            x[k] + lx[k] <= x[i] + M * (1 - right[i, k]),
            name=f"right_sep_{i}_{k}"
        )

        model.addConstr(
            y[i] + wy[i] <= y[k] + M * (1 - front[i, k]),
            name=f"front_sep_{i}_{k}"
        )
        model.addConstr(
            y[k] + wy[k] <= y[i] + M * (1 - back[i, k]),
            name=f"back_sep_{i}_{k}"
        )

        model.addConstr(
            z[i] + hz[i] <= z[k] + M * (1 - below[i, k]),
            name=f"below_sep_{i}_{k}"
        )
        model.addConstr(
            z[k] + hz[k] <= z[i] + M * (1 - above[i, k]),
            name=f"above_sep_{i}_{k}"
        )

        model.addConstr(left[i, k] + right[i, k] <= 1, name=f"x_mutex_{i}_{k}")
        model.addConstr(front[i, k] + back[i, k] <= 1, name=f"y_mutex_{i}_{k}")
        model.addConstr(below[i, k] + above[i, k] <= 1, name=f"z_mutex_{i}_{k}")

    model.setObjective(0.0, GRB.MINIMIZE)

    aux = {
        "base_fitness": base_fitness,
        "trying_customer": trying_customer,
        "saved_outsource_cost": saved_outsource_cost,
        "extra_fuel_cost": extra_fuel_cost,
        "estimated_rescued_fitness_if_success": base_fitness - saved_outsource_cost + extra_fuel_cost,

        "old_route": old_route,
        "local_rescued_route": local_rescued_route,
        "full_planned_route": full_planned_route,
        "remaining_after_fail": remaining_after_fail,
        "fail_index_in_planned_route": fail_index,

        "old_dist": old_dist,
        "local_rescued_dist": local_rescued_dist,
        "full_planned_dist": full_planned_dist,

        "vi": vi,
        "outsource_fee": outsource_fee,
    }

    var_dict = {
        "x": x,
        "y": y,
        "z": z,
        "placed": placed,
        "o": o,
        "lx": lx,
        "wy": wy,
        "hz": hz,
        "rescue_success": rescue_success,
        "left": left,
        "right": right,
        "front": front,
        "back": back,
        "below": below,
        "above": above,
    }

    return model, var_dict, aux
# =========================================================
# main
# =========================================================
if __name__ == "__main__":
    case_data = read_single_truck_case(base_dir, case_file)
    coords = read_generated_coords(base_dir, coord_file)

    model, var_dict, aux = build_single_truck_packing_model(case_data, coords, C1=5.41)

    print("case_id =", case_data["case_id"])
    print("trying_customer =", case_data["trying_customer"])
    print("old_route =", aux["old_route"])
    print("local_rescued_route =", aux["local_rescued_route"])
    print("full_planned_route =", aux["full_planned_route"])
    print("remaining_after_fail =", aux["remaining_after_fail"])
    print("fail_index_in_planned_route =", aux["fail_index_in_planned_route"])
    print("old_dist =", aux["old_dist"])
    print("local_rescued_dist =", aux["local_rescued_dist"])
    print("saved_outsource_cost =", aux["saved_outsource_cost"])
    print("extra_fuel_cost =", aux["extra_fuel_cost"])
    print("estimated_rescued_fitness_if_success =", aux["estimated_rescued_fitness_if_success"])

    model.setParam("TimeLimit", 30)
    model.optimize()

    print("model status =", model.Status)

    if model.Status == GRB.OPTIMAL:
        print("Packing feasible: YES")
        print("Rescue success: YES")
        print("rescued fitness =", aux["estimated_rescued_fitness_if_success"])

        for i, item in enumerate(case_data["items"]):
            chosen_ori = None
            for r in range(6):
                if var_dict["o"][i, r].X > 0.5:
                    chosen_ori = r
                    break

            print(
                f"item {i} | "
                f"cust={item['customer_id']} cargo={item['cargo_id']} | "
                f"x={var_dict['x'][i].X:.3f}, "
                f"y={var_dict['y'][i].X:.3f}, "
                f"z={var_dict['z'][i].X:.3f}, "
                f"lx={var_dict['lx'][i].X:.3f}, "
                f"wy={var_dict['wy'][i].X:.3f}, "
                f"hz={var_dict['hz'][i].X:.3f}, "
                f"ori={chosen_ori}"
            )

    elif model.Status == GRB.INFEASIBLE:
        print("Packing feasible: NO")
        print("Rescue success: NO")
        print("this case cannot be packed under current model")

    elif model.Status == GRB.TIME_LIMIT:
        print("Time limit reached")
        print("best status not proven optimal")
        if model.SolCount > 0:
            print("But at least one feasible solution was found")
        else:
            print("No feasible solution found within time limit")

    else:
        print("Packing result unclear")
        print("status code =", model.Status)