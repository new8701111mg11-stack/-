# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from math import sqrt
import random
import math


# 寫入 CSV 檔
def write_to_csv(goods_data, service_areas, num_customers, routes, arcs):

    # 客戶貨物數量
    customer_goods_stats = []

    for customer_id in range(1, num_customers + 1):
        # 篩選該客戶的所有貨物資料
        customer_rows = [row for row in goods_data if row[0] == customer_id]
        
        count = len(customer_rows)
        total_volume = sum(row[2] for row in customer_rows)  # row[2] 是 volume
        
        customer_goods_stats.append([customer_id, count, total_volume])

    # 建立 DataFrame 並輸出
    count_df = pd.DataFrame(customer_goods_stats, columns=['客戶', '貨物數', '總才積'])
    count_df.to_csv("customerInfo.csv", index=False, encoding='utf-8-sig')

    # 貨物明細
    goods_df = pd.DataFrame(goods_data, columns=["客戶", "貨物", "材積", "長", "寬", "高", "方向1", "方向2", "方向3", "方向4", "方向5", "方向6", "脆弱性"])
    goods_df.to_csv("goods.csv", index=False, encoding='utf-8-sig')

    # 服務區域
    service_area_columns = ["客戶"] + [f"服務區域{i+1}" for i in range(len(service_areas[0]) - 1)]
    service_df = pd.DataFrame(service_areas, columns=service_area_columns)
    service_df.to_csv("serviceArea.csv", index=False, encoding='utf-8-sig')

    #車輛路線
    max_len = max(len(route) for route in routes.values())
    route_rows = []
    for area_id, route in routes.items():
        padded_route = route + [""] * (max_len - len(route))  # 補空白格
        route_rows.append([area_id] + padded_route)

    route_columns = ["區域"] + [f"節點{i+1}" for i in range(max_len)]
    routes_df = pd.DataFrame(route_rows, columns=route_columns)
    routes_df.to_csv("routes.csv", index=False, encoding='utf-8-sig')

    #arc
    arc_df = pd.DataFrame(arcs, columns=["區域", "起點", "終點"])
    arc_df.to_csv("routeArcs.csv", index=False, encoding="utf-8-sig")

# 生成貨物資訊
item_dimensions = [
    [57, 51, 54], [92, 81, 35], [81, 33, 28], [98, 122, 46], [92, 12, 56],
    [83, 100, 30], [73, 57, 60], [91, 71, 32], [40, 39, 60],
    [82, 67, 30], [33, 25, 23], [86, 38, 22], [84, 112, 43], [47, 36, 45], [81, 78, 47], [50, 51, 48]
]
def generate_customer_goods(num_customers, max_goods_per_customer = 3):
    data = []
    for customer_id in range(1, num_customers + 1):
        num_goods = random.randint(1, max_goods_per_customer)
        for good_id in range(1, num_goods + 1):
            length, width, height = random.choice(item_dimensions)
            volume = length * width * height
            orientation_flags = [0] * 6
            first_one = random.randint(0, 5)
            orientation_flags[first_one] = 1
            for i in range(6):
                if i != first_one:
                    orientation_flags[i] = random.randint(0, 1)
            fragile = 1 if random.random() < 0.2 else 0
            data.append([customer_id, good_id, volume, length, width, height] + orientation_flags + [fragile])
    return data

def generate_service_areas(num_customers, num_areas):
    # 初始化顧客分配
    customer_ids = list(range(1, num_customers + 1))
    random.shuffle(customer_ids)

    base_count = num_customers // num_areas
    remainder = num_customers % num_areas

    service_areas = []
    area_to_customers = {}  # 儲存每區有哪些顧客

    start = 0
    for area_id in range(num_areas):
        count = base_count + (1 if area_id < remainder else 0)
        assigned = customer_ids[start:start + count]
        area_to_customers[area_id] = assigned
        for cust_id in assigned:
            zone_flags = [0] * num_areas
            zone_flags[area_id] = 1
            service_areas.append([cust_id] + zone_flags)
        start += count

    # 依區域選出 1/3 顧客進行鄰居擴張
    for area_id in range(num_areas):
        customers = area_to_customers[area_id]
        num_overlap = max(1, math.ceil(len(customers) / 3))
        selected = random.sample(customers, num_overlap)

        for cust_id in selected:
            row = next(r for r in service_areas if r[0] == cust_id)
            # 加入右邊區域（除了最後一區）
            if area_id < num_areas - 1:
                row[1 + area_id + 1] = 1
            # 最後一區，加左邊區域
            elif area_id == num_areas - 1 and area_id - 1 >= 0:
                row[1 + area_id - 1] = 1

    # 最後依顧客 ID 排序
    service_areas.sort(key=lambda r: r[0])
    return service_areas

def generate_vehicle_routes(service_areas, num_areas, shuffle=True):
    routes = {}
    for area_id in range(num_areas):
        route = [0]  # 起始點
        customers_in_area = []
        
        for row in service_areas:
            cust_id = row[0]
            zone_flags = row[1:]
            if zone_flags[area_id] == 1:
                customers_in_area.append(cust_id)

        if shuffle:
            random.shuffle(customers_in_area)

        route.extend(customers_in_area)
        route.append(0) #終點
        routes[area_id] = route

    return routes

# def generate_arcs(routes):
#     arc_list = []
#     for area_id, route in routes.items():
#         for i in range(len(route) - 1):
#             from_node = route[i]
#             to_node = route[i + 1]
#             arc_list.append([area_id, from_node, to_node])

#     return arc_list
def generate_arcs(routes):
    arc_list = []
    for area_id, route in routes.items():
        for i in range(len(route)):
            for j in range(i + 1, len(route)):
                from_node = route[i]
                to_node = route[j]
                if from_node != to_node:
                    arc_list.append([area_id, from_node, to_node])
    return arc_list
# main function
num_customers = 60
num_areas = 4
goods_data = generate_customer_goods(num_customers)
service_areas = generate_service_areas(num_customers, num_areas)
routes = generate_vehicle_routes(service_areas, num_areas=4)
arcs = generate_arcs(routes)
write_to_csv(goods_data, service_areas, num_customers, routes, arcs)
print("所有資料已輸出為 CSV 檔案 ✅")