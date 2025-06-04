# import pickle
# import numpy as np
# from collections.abc import Callable


# def load_table(name: str) -> np.ndarray:
#     with open(f"tables/{name}.pkl", "rb") as file:
#         return pickle.load(file)


# def save_table(name: str, table: np.ndarray):
#     with open(f"tables/{name}.pkl", "wb") as file:
#         pickle.dump(table, file)


# def get_tables(tables_names: list, tables_kwargs: list, generate_table_func: Callable) -> list:
#     tables = []
#     for phase, tables_names in enumerate(tables_names):
#         phase_tables = []
#         for idx, name in enumerate(tables_names):
#             try:
#                 table = load_table(name)
#             except FileNotFoundError:
#                 table = generate_table_func(phase, idx, **tables_kwargs[phase][idx])
#                 save_table(name, table)
#             phase_tables.append(table)
#         tables.append(phase_tables)
#     return tables
