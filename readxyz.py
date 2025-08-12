import numpy as np

def read_xyz_file(file_path):
    """读取XYZ文件，返回仅含C原子的坐标和元素符号列表"""
    with open(file_path, 'r') as f:
        lines = f.readlines()
        num_atoms = int(lines[0].strip())
        atoms = []
        for line in lines[2: 2+num_atoms]:
            parts = line.strip().split()
            element = parts[0]
            if element == 'C':  # 仅保留C原子
                x, y, z = list(map(float, parts[1:4]))
                atoms.append((element, x, y, z))
    # 提取C原子坐标和元素列表
    elements = [atom[0] for atom in atoms]  # 全部为'C'
    positions = [atom[1:] for atom in atoms]
    return positions, elements

def calculate_distances(positions):
    """计算所有原子对之间的欧氏距离"""
    n = len(positions)
    distance_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            dx = positions[i][0] - positions[j][0]
            dy = positions[i][1] - positions[j][1]
            dz = positions[i][2] - positions[j][2]
            distance = np.sqrt(dx**2 + dy**2 + dz**2)
            distance_matrix[i][j] = distance
            distance_matrix[j][i] = distance
    return distance_matrix

def generate_adjacency_matrix(distance_matrix, threshold=1.5):  #threshold参数可能需要修改
    n = len(distance_matrix)
    adjacency_matrix = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if distance_matrix[i][j] < threshold:  #检查距离
                adjacency_matrix[i][j] = 1
                adjacency_matrix[j][i] = 1
    return adjacency_matrix

if __name__ == "__main__":
    file_path = "molecule.xyz"
    positions, elements = read_xyz_file(file_path)  # 此处的elements全部为C
    distance_matrix = calculate_distances(positions)
    adjacency_matrix = generate_adjacency_matrix(distance_matrix)  # 维度为C原子数×C原子数
    
    print("邻接矩阵（仅C-C键，维度为{}×{}）:".format(len(positions), len(positions)))
    print(adjacency_matrix)



