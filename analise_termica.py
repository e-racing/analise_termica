from math import sqrt, exp, log
import numpy as np
import pandas as pd
import plotly.express as px

class Segmento:
    def __init__(self, diametro, espaçamento_transversal, espaçamento_longitudinal, l, v_in, numero_de_celulas_longitudinal, 
                 numero_de_celulas_transversal, temperatura_superficie):
        """
        diametro -> diametro da célula
        st -> passo transversal entre células
        sl -> passo longitudinal entre células
        l -> comprimento da célula
        v_in -> velocidade de entrada do ar
        """
        self.diametro = diametro
        self.st = espaçamento_transversal + diametro
        self.sl = espaçamento_longitudinal + diametro
        self.l = l
        self.v_in = v_in
        # Pre-calculate sd for reuse
        self.sd = sqrt((self.st / 2) ** 2 + self.sl ** 2)
        self.numero_de_celulas_longitudinal = numero_de_celulas_longitudinal
        self.numero_de_celulas_transversal = numero_de_celulas_transversal
        self.area_transferencia_calor = self.numero_de_celulas_longitudinal * self.numero_de_celulas_transversal * self.diametro * self.l
        self.temperatura_superficie = temperatura_superficie


    def VelocidadeMaxima(self):
        # Calculate areas of flow
        at = (self.st - self.diametro) * self.l
        ad = (self.sd - self.diametro) * self.l

        # Determine maximum velocity based on area conditions
        if 2 * ad > at:
            v_max = self.st * self.v_in / (self.st - self.diametro)
        else:
            v_max = self.st * self.v_in / (2 * (self.sd - self.diametro))

        return v_max

    def NumeroDeReynolds(self, densidade, viscosidade):
        """
        Calculates Reynolds number using maximum velocity.
        """
        v_max = self.VelocidadeMaxima()  # Only call once
        re = (densidade * v_max * self.diametro) / viscosidade
        return re

    def NumeroDeNusselt(self, densidade, viscosidade, pr_superficie, pr_medio):
        """
        Calculates Nusselt number using Reynolds number and corrects for the number of rows.
        """
        re = self.NumeroDeReynolds(densidade, viscosidade)
        
        if re <= 500:
            nu = 1.04 * re ** 0.4 * pr_medio ** 0.36 * (pr_medio / pr_superficie) ** 0.25
        elif re <= 1000:
            nu = 0.71 * re ** 0.5 * pr_medio ** 0.36 * (pr_medio / pr_superficie) ** 0.25
        elif re <= 2e5:
            nu = 0.35 * (self.st / self.sl) ** 0.2 * re ** 0.6 * pr_medio ** 0.36 * (pr_medio / pr_superficie) ** 0.25
        elif re < 2e6:
            nu = 0.031 * (self.st / self.sl) ** 0.2 * re ** 0.8 * pr_medio ** 0.36 * (pr_medio / pr_superficie) ** 0.25
        else:
            nu = 0

        # Known values for number of rows and correction factors
        fileiras_conhecidas = [1, 2, 3, 4, 5, 7, 10, 13]
        fatores_conhecidos = [0.64, 0.76, 0.84, 0.89, 0.93, 0.96, 0.98, 0.99]

        # Check if rows need correction, using interpolation if necessary
        if self.numero_de_celulas_transversal in fileiras_conhecidas:
            indice = fileiras_conhecidas.index(self.numero_de_celulas_transversal)
            fator_corrigido = fatores_conhecidos[indice]
        else:
            fator_corrigido = np.interp(self.numero_de_celulas_transversal, fileiras_conhecidas, fatores_conhecidos)

        nu_corrigido = fator_corrigido * nu

        return nu, nu_corrigido

    def CoeficienteDeTransferenciaDeCalor(self, densidade, viscosidade, pr_medio, pr_superficie, condutividade_termica):
        """
        Calculates the heat transfer coefficient.
        """
        _, nu_corrigido = self.NumeroDeNusselt(densidade, viscosidade, pr_medio, pr_superficie)
        coeficiente_h = nu_corrigido * (condutividade_termica / self.diametro)
        return coeficiente_h

    def VazaoMassica(self, densidade):
        """
        Calculates mass flow rate.
        """
        return densidade * self.v_in * self.numero_de_celulas_transversal * self.st * self.l
    
    def TemperaturaSaidaFluido(self, densidade, temperatura_superficie, temperatura_entrada, coeficiente_h):

        vazao_massica = self.VazaoMassica(densidade)

        temp_saida = temperatura_superficie - (temperatura_superficie - temperatura_entrada) * exp(-self.area_transferencia_calor * coeficiente_h / (vazao_massica * coeficiente_h))

        return temp_saida
    
    def DiferencaMediaLogaritmicaDeTemperatura(self, temperatura_superficie, temperatura_saida, temperatura_entrada):

        dif_superficie_saida = temperatura_superficie - temperatura_saida

        dif_superficie_entrada = temperatura_superficie - temperatura_entrada

        diferenca_media_log_temp = (dif_superficie_saida - dif_superficie_entrada) / log(dif_superficie_saida / dif_superficie_entrada)

        return diferenca_media_log_temp
    
    def TransferenciaDeCalor(self, coeficiente_h, diferenca_media_log_temp):

        transferencia_de_calor = coeficiente_h * self.area_transferencia_calor * diferenca_media_log_temp

        return transferencia_de_calor 

    def QuedaDePressao(self, coeficiente_de_atrito, fator_de_correcao, densidade, v_max):
        
        queda_de_pressao = self.numero_de_celulas_longitudinal * coeficiente_de_atrito * fator_de_correcao * densidade * pow(v_max, 2) / 2

        return queda_de_pressao

    def PotenciaNecessaria(self, vazao_massica, queda_pressao, densidade): # FONTE ?

        potencia = vazao_massica * queda_pressao / densidade

        return potencia

diametro_celula = 0.018
espacamento_transversal = 0.001
espacamento_longitudinal = 0.001
comprimento_celula = 0.065
velocidade_entrada_ar = 5
numero_celulas_longitudinal = 5
numero_celulas_transversal = 6
temperatura_superficie = 40
densidade = 1.145
viscosidade = 1.895e-5
pr_superficie = 0.7255
pr_medio = 0.7268
condutividade_termica = 0.02625
temperatura_entrada = 35
coeficiente_de_atrito = 0.16
fator_correcao = 1

segmento = Segmento(diametro_celula, espacamento_transversal, espacamento_longitudinal, comprimento_celula, velocidade_entrada_ar, 
                    numero_celulas_longitudinal, numero_celulas_transversal, temperatura_superficie)


vel_max = segmento.VelocidadeMaxima()
numero_reynolds = segmento.NumeroDeReynolds(densidade, viscosidade)
numero_nusselt, numero_nusselt_corrigido = segmento.NumeroDeNusselt(densidade, viscosidade, pr_superficie, pr_medio)
coeficiente_h = segmento.CoeficienteDeTransferenciaDeCalor(densidade, viscosidade, pr_medio, pr_superficie, condutividade_termica)
vazao_massica = segmento.VazaoMassica(densidade)
temp_saida_fluido = segmento.TemperaturaSaidaFluido(densidade, temperatura_superficie, temperatura_entrada, coeficiente_h)
dif_med_log_temp = segmento.DiferencaMediaLogaritmicaDeTemperatura(temperatura_superficie, temp_saida_fluido, temperatura_entrada)
transf_calor = segmento.TransferenciaDeCalor(coeficiente_h, dif_med_log_temp)
queda_pressao = segmento.QuedaDePressao(coeficiente_de_atrito, fator_correcao, densidade, vel_max)
potencia_necessaira = segmento.PotenciaNecessaria(vazao_massica, queda_pressao, densidade)

print(f"""
      vel_max: {vel_max} m/s
      reynolds: {numero_reynolds} 
      nusselt: {numero_nusselt}
      numero_nusselt_corrigido: {numero_nusselt_corrigido} 
      coeficiente_h: {coeficiente_h} W / m^2 K
      vazao massica: {vazao_massica} Kg/s
      temperatura de saida do fluido: {temp_saida_fluido} C
      diferenca media logartmica de temperatura: {dif_med_log_temp} C
      transferencia de calor: {transf_calor} W
      queda de pressao: {queda_pressao} Pa
      potencia necessaria: {potencia_necessaira} W
    
""")



# seg = Segmento(0.018, 0.001, 0.001, 0.065, 5)
# v_max = seg.VelocidadeMaxima()
# re = seg.NumeroDeReynolds(densidade, viscosidade)
# nu, nu_corrigido = seg.NumeroDeNusselt(densidade, viscosidade, pr_medio, pr_superficie, self.numero_de_celulas_transversal)
# coeficiente_h = seg.CoeficienteDeTransferenciaDeCalor(densidade, viscosidade, pr_medio, pr_superficie, self.numero_de_celulas_transversal, condutividade_termica)
# vazao_massica = seg.VazaoMassica(densidade, numero_de_celulas_transversal)
# temp_saida, diferenca_media_log_temp, transf_de_calor = seg.TransferenciaDeCalor(densidade, numero_de_celulas_transversal, self.numero_de_celulas_transversal, temperatura_superficie, temp_entr, viscosidade, pr_medio, pr_superficie, condutividade_termica, calor_especifico)

# espacamento_transversal = np.arange(0.001, 0.021, 0.001)
# espacamento_longitudinal = np.arange(0.001, 0.021, 0.001)
# transf_de_cal = []

# for espacamento_tr in espacamento_transversal:
#   for espacamento_long in espacamento_longitudinal:

#     seg = Segmento(0.018, espacamento_tr, espacamento_long, 0.065, 5)
#     v_max = seg.VelocidadeMaxima()
#     re = seg.NumeroDeReynolds(densidade, viscosidade)
#     nu, nu_corrigido = seg.NumeroDeNusselt(densidade, viscosidade, pr_medio, pr_superficie, self.numero_de_celulas_transversal)
#     coeficiente_h = seg.CoeficienteDeTransferenciaDeCalor(densidade, viscosidade, pr_medio, pr_superficie, self.numero_de_celulas_transversal, condutividade_termica)
#     vazao_massica = seg.VazaoMassica(densidade, numero_de_celulas_transversal)
#     temp_saida, diferenca_media_log_temp, transf_de_calor = seg.TransferenciaDeCalor(densidade, numero_de_celulas_transversal, self.numero_de_celulas_transversal, temperatura_superficie, temp_entr, viscosidade, pr_medio, pr_superficie, condutividade_termica, calor_especifico)
#     transf_de_cal.append(transf_de_calor)


# data = {
#     'Espacamento Transversal': np.repeat(espacamento_transversal, len(espacamento_longitudinal)),
#     'Espacamento Longitudinal': np.tile(espacamento_longitudinal, len(espacamento_transversal)),
#     'Transf de Cal': transf_de_cal
# }

# df = pd.DataFrame(data)

# # Criando o gráfico de dispersão
# fig = px.scatter(df, x='Espacamento Transversal', y='Espacamento Longitudinal', color='Transf de Cal',
#                  color_continuous_scale='Viridis', title='Gráfico de Dispersão de Transferência de Calor')

# fig.show()
