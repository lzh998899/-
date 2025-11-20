import numpy as np
from scipy.integrate import quad

# Henyey-Greenstein相位函数
def p(cos_theta, g):
    return (1 - g**2) / (4 * np.pi * (1 + g**2 - 2 * g * cos_theta)**(3/2))

# RTE辐射传输方程的计算
def RTE_solver(I0, mu_a, mu_s, d, g):
    # 对散射角度积分，计算散射对强度的贡献
    scattering_integral, _ = quad(lambda cos_theta: p(cos_theta, g), -1, 1)
    
    # 考虑吸收和散射后的最终强度计算  mu-s 后面跟着的
    I_d = I0 * np.exp(-(mu_a + mu_s) * d) + (mu_s/ (4 * np.pi) ) * scattering_integral * I0 * d
    return I_d

# 组织参数列表
organs = [
    {
        'name': 'cornea',
        'I0': 3.144654,  # mW/cm^2, 入射辐照度
        'mu_a': 0.,  # mm^-1, 吸收系数
        'mu_s': 0.0016,  # mm^-1, 散射系数
        'd': 0.25,       # mm, 角膜厚度
        'g': 0.9         # Henyey-Greenstein相位函数的各向异性因子
    },
    {
        'name': 'aqueous',
        'mu_a': 0.0078,  # mm^-1, 吸收系数
        'mu_s': 0.0060,  # mm^-1, 散射系数
        'd': 1.60023,    # mm, 房水厚度
        'g': 0.9         # Henyey-Greenstein相位函数的各向异性因子
    },
    {
        'name': 'lens',
        'mu_a': 0.0005,  # mm^-1, 吸收系数
        'mu_s': 0.0023,  # mm^-1, 散射系数
        'd': 2.77023,    # mm, 晶体厚度
        'g': 0.9         # Henyey-Greenstein相位函数的各向异性因子
    },
    {
        'name': 'vitreous',
        'mu_a': 0.048,   # mm^-1, 吸收系数
        'mu_s': 0.0003,  # mm^-1, 散射系数
        'd': 3.17659,    # mm, 玻璃体厚度
        'g': 0.9         # Henyey-Greenstein相位函数的各向异性因子
    },
    {
        'name': 'retina',
        'mu_a': 0.02,    # mm^-1, 吸收系数
        'mu_s': 5,       # mm^-1, 散射系数
        'd': 0.15,       # mm, 玻璃体厚度
        'g': 0.9         # Henyey-Greenstein相位函数的各向异性因子
    },
    {
        'name': 'choroid',
        'mu_a': 0.1,     # mm^-1, 吸收系数
        'mu_s': 20,      # mm^-1, 散射系数
        'd': 0.1,        # mm, 玻璃体厚度
        'g': 0.9         # Henyey-Greenstein相位函数的各向异性因子
    }
]

# 初始化I0为角膜的入射辐照度
current_I0 = organs[0]['I0']

# 输出初始入射辐照度
print(f"初始入射辐照度（角膜）：{current_I0} mW/cm²")

# 遍历每个组织，计算经过后的衰减辐照度
for organ in organs:
    # 计算经过当前组织后的衰减辐照度
    attenuated_intensity = RTE_solver(current_I0, organ['mu_a'], organ['mu_s'], organ['d'], organ['g'])
    
    # 输出经过当前组织后的衰减辐照度
    print(f"经过{organ['name']}后的衰减辐照度为: {attenuated_intensity:.4f} mW/cm²")
    
    # 将当前组织的衰减辐照度作为下一个组织的入射辐照度
    current_I0 = attenuated_intensity
