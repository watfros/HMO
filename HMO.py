import numpy as np

class HMO_Calculator:
    def __init__(self, alpha=0.0, beta=-1.0):
        self.alpha = alpha  # 库仑积分参数α
        self.beta = beta    # 共振积分参数β
        self.bond_orders = None
        self.charge_density = None
        self.free_valence = None
        self.molecular_orbitals = None
        self.orbital_energies = None
        self.delocalization_energy = None  # 离域能

    def build_huckel_matrix(self, adjacency_matrix):
        """构建Hückel矩阵"""
        n = len(adjacency_matrix)
        H = np.zeros((n,n))
        for i in range(n):
            H[i,i] = self.alpha  # 对角元素设为α
            for j in range(n):
                if adjacency_matrix[i][j] == 1:
                    H[i,j] = self.beta  # 相邻原子设为β
        return H

    def solve_mo(self, H_matrix, num_electrons):
        """解分子轨道"""
        eigenvalues, eigenvectors = np.linalg.eigh(H_matrix)
        sorted_idx = np.argsort(eigenvalues)
        self.orbital_energies = eigenvalues[sorted_idx]
        self.molecular_orbitals = eigenvectors[:, sorted_idx]
        self.occupied = num_electrons // 2  # 占据轨道数
        
    def calculate_properties(self):
        """计算量子化学指数"""
        n = self.molecular_orbitals.shape[0]
        # 电荷密度计算
        self.charge_density = np.zeros(n)
        for atom in range(n):
            density = 0.0
            for mo in range(self.occupied):
                density += 2 * (self.molecular_orbitals[atom, mo])**2
            self.charge_density[atom] = density
        
        # 键级计算  
        self.bond_orders = np.zeros((n,n))
        for i in range(n):
            for j in range(i+1, n):
                p = 0.0
                if adjacency[i][j] == 1:  # 仅累加相邻原子键级[4][7]
                    for mo in range(self.occupied):
                        p += 2 * self.molecular_orbitals[i,mo] * self.molecular_orbitals[j,mo]
                    self.bond_orders[i,j] = p
                    self.bond_orders[j,i] = p  # 确保矩阵对称
        
        # 自由价计算：根据邻接矩阵筛选相邻原子键级
        self.free_valence = np.zeros(n)
        for atom in range(n):
            bond_sum = 0.0
            for neighbor in range(n):
                if adjacency[atom][neighbor] == 1:  # 仅累加相邻原子键级[4][7]
                    bond_sum += self.bond_orders[atom, neighbor]
            self.free_valence[atom] = np.sqrt(3) - bond_sum  # 碳原子N_max=√3[4]

        # 离域能计算
        # 计算离域体系的总能量
        total_delocalized_energy = 2 * np.sum(self.orbital_energies[:self.occupied])
        #print(total_delocalized_energy)
        # 计算定域体系的总能量（假设每个双键单独存在）
        # 双键的数量等于占据数
        num_double_bonds = self.occupied

        # 定域体系的总能量（每个双键能量为 alpha + beta）
        total_localized_energy = 2 * num_double_bonds * (self.alpha + self.beta)
        print(total_localized_energy)
        # 离域能为两者的差值
        self.delocalization_energy = total_delocalized_energy - total_localized_energy

    def print_results(self):
        """输出结果"""
        print("分子轨道能量：")
        for i, energy in enumerate(self.orbital_energies):
            print(f"轨道 {i+1}: {energy:.4f}")
        
        print("\n分子轨道组成：")
        for i in range(len(self.molecular_orbitals)):
            print(f"轨道 {i+1} 系数：{self.molecular_orbitals[:,i].round(4)}")
        
        print("\n电荷密度：")
        for i, rho in enumerate(self.charge_density):
            print(f"原子 {i+1}: {rho:.4f}")
        
        print("\n键级：")
        n = len(self.bond_orders)
        for i in range(n):
            for j in range(i+1, n):
                if self.bond_orders[i,j] > 1e-4:
                    print(f"原子 {i+1}-{j+1}: {self.bond_orders[i,j]:.4f}")
        
        print("\n自由价：")
        for i, f in enumerate(self.free_valence):
            print(f"原子 {i+1}: {f:.4f}")

        print(f"\n离域能: {self.delocalization_energy:.4f}")

# 示例使用：分子
if __name__ == "__main__":
    # 邻接矩阵
    adjacency = [
      [0, 1, 0, 0],
      [1, 0, 1, 0],
      [0, 1, 0, 1],
      [0, 0, 1, 0]
    ]
    
    hmo = HMO_Calculator(alpha=0, beta=-1)
    H = hmo.build_huckel_matrix(adjacency)
    # 修改π电子数目
    hmo.solve_mo(H, num_electrons=4)  
    hmo.calculate_properties()
    hmo.print_results()