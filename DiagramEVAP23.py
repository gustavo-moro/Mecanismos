import numpy as np
import matplotlib.pyplot as plt

# Altura máxima
h = 10
Rb = 17.5

# Definir o intervalo de theta (em graus)
theta_degrees = np.linspace(0, 360, 1000)

# Inicializar os vetores para os resultados
E = np.zeros_like(theta_degrees)
V = np.zeros_like(theta_degrees)
A = np.zeros_like(theta_degrees)

# Parâmetros do movimento
subida = 55  # Ângulo de subida em graus
espera = 10  # Ângulo de espera superior em graus
descida = 55  # Ângulo de descida em graus
ciclo = subida + espera + descida  # Ciclo completo

# Função perfil 2-3 para subida e descida
def perfil_23(theta, theta0, beta_degrees, h0, h_max):
    """
    Calcula o perfil de posição, velocidade e aceleração usando o perfil 2-3.
    """
    # Converter beta de graus para radianos
    beta = np.radians(beta_degrees)
    delta_h = h_max - h0  # Diferença de altura (pode ser negativa para descida)
    t = (np.radians(theta) - np.radians(theta0)) / beta  # Normalizar o ângulo
    E = h0 + delta_h * (3 * t**2 - 2 * t**3)
    V = (delta_h / beta) * (6 * t - 6 * t**2)  # Velocidade
    A = (delta_h / beta**2) * (6 - 12 * t)  # Aceleração
    return E, V, A

# Preencher as curvas para os três ciclos
for i in range(3):  # Três ciclos completos
    # Definir os limites de cada etapa
    theta0_subida = i * ciclo
    theta0_espera = theta0_subida + subida
    theta0_descida = theta0_espera + espera

    # Subida
    mask_subida = (theta_degrees >= theta0_subida) & (theta_degrees < theta0_subida + subida)
    E[mask_subida], V[mask_subida], A[mask_subida] = perfil_23(
        theta_degrees[mask_subida], theta0_subida, subida, 0, h
    )

    # Espera superior (manter altura máxima)
    mask_espera = (theta_degrees >= theta0_espera) & (theta_degrees < theta0_espera + espera)
    E[mask_espera] = h
    V[mask_espera] = 0
    A[mask_espera] = 0

    # Descida
    mask_descida = (theta_degrees >= theta0_descida) & (theta_degrees < theta0_descida + descida)
    E[mask_descida], V[mask_descida], A[mask_descida] = perfil_23(
        theta_degrees[mask_descida], theta0_descida, descida, h, 0
    )

# Calcular o ângulo de pressão para todos os pontos
with np.errstate(divide='ignore', invalid='ignore'):  # Ignorar divisões por zero
    phi = np.degrees(np.arctan(V / (E + np.sqrt(Rb**2))))

# Dados do problema
rpm = 120  # Velocidade de rotação em rpm

# Converter rpm para rad/s
omega = 2 * np.pi * rpm / 60  # Velocidade angular em rad/s

# Encontrar os índices dos máximos para cada variável
idx_max_v = np.argmax(V)  # Índice do valor máximo de velocidade
idx_max_a = np.argmax(A)  # Índice do valor máximo de aceleração
idx_max_phi = np.argmax(phi)  # Índice do valor máximo de ângulo de pressão

# Obter os valores máximos
max_v = V[idx_max_v]
max_v_seguidor = max_v * omega
max_a = A[idx_max_a]
max_a_seguidor = max_a * omega
max_phi = phi[idx_max_phi]

# Obter os ângulos correspondentes aos máximos
theta_max_v = theta_degrees[idx_max_v]
theta_max_a = theta_degrees[idx_max_a]
theta_max_phi = theta_degrees[idx_max_phi]

# Exibir os resultados
print(f"Máximo de Velocidade: {max_v:.2f} mm/rad em θ = {theta_max_v:.2f}°")
print(f"Velocidade máxima do seguidor: {max_v_seguidor:.2f} mm/s em θ = {theta_max_v:.2f}°")
print(f"Máximo de Aceleração: {max_a:.2f} mm/rad² em θ = {theta_max_a:.2f}°")
print(f"Aceleração máxima do seguidor: {max_a_seguidor:.2f} mm/rad² em θ = {theta_max_a:.2f}°")
print(f"Máximo de Ângulo de Pressão: {max_phi:.2f}° em θ = {theta_max_phi:.2f}°")



# Plotar os gráficos
fig, axs = plt.subplots(3, 1, figsize=(10, 10))

# Ponto de referência para marcação no eixo X
ticks_personalizados = [55, 65, 120, 175, 185, 240, 295, 305, 360]

axs[0].plot(theta_degrees, E, label="Posição E(θ)")
axs[0].set_ylabel("Posição E(θ) (mm)")
axs[0].grid()
axs[0].legend()

axs[1].plot(theta_degrees, V, label="Velocidade V(θ)", color="orange")
axs[1].set_ylabel("Velocidade V(θ) (mm/rad)")
axs[1].grid()
axs[1].legend()

axs[2].plot(theta_degrees, A, label="Aceleração A(θ)", color="green")
axs[2].set_xlabel("Ângulo θ (graus)")
axs[2].set_ylabel("Aceleração A(θ) (mm/rad²)")
axs[2].grid()
axs[2].legend()

# Definir o formato de coordenadas para cada subgráfico
def mostrar_coordenadas(x, y):
    return f"x={x:.2f}, y={y:.2f}"  # Formato com duas casas decimais

# Aplicar a função em cada eixo
for ax in axs:
    ax.format_coord = mostrar_coordenadas

# Ajustar os ticks personalizados para todos os subgráficos
for ax in axs:
    ax.set_xticks(ticks_personalizados)
    ax.set_xticklabels(ticks_personalizados)  # Garantir que os labels reflitam os valores
    ax.set_xlim([0, 360])  # Limitar o eixo X de 0 a 360 graus

# Ajustar o layout e exibir
plt.tight_layout()
plt.show()

# Criar gráfico separado para o Ângulo de Pressão
fig2, ax2 = plt.subplots(figsize=(10, 4))

ax2.plot(theta_degrees, phi, label="Ângulo de Pressão ϕ(θ)", color="red")
ax2.set_xlabel("Ângulo θ (graus)")
ax2.set_ylabel("Ângulo de Pressão ϕ(θ) (graus)")
ax2.grid()
ax2.legend()

# Aplicar a função em cada eixo
ax2.format_coord = mostrar_coordenadas

# Ajustar os ticks personalizados para todos os subgráficos
ax2.set_xticks(ticks_personalizados)
ax2.set_xticklabels(ticks_personalizados)  # Garantir que os labels reflitam os valores
ax2.set_xlim([0, 360])  # Limitar o eixo X de 0 a 360 graus

# Exibir o gráfico do ângulo de pressão
plt.tight_layout()
plt.show()
