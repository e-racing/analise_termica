from math import sqrt, exp, log
import numpy as np
import pandas as pd
import plotly.express as px

class Segmento:
    def __init__(self, diametro, espaçamento_transversal, espaçamento_longitudinal, l, v_in):
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

    def NumeroDeNusselt(self, densidade, viscosidade, pr_superficie, pr_medio, numero_fileiras):
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
        if numero_fileiras in fileiras_conhecidas:
            indice = fileiras_conhecidas.index(numero_fileiras)
            fator_corrigido = fatores_conhecidos[indice]
        else:
            fator_corrigido = np.interp(numero_fileiras, fileiras_conhecidas, fatores_conhecidos)

        nu_corrigido = fator_corrigido * nu

        return nu, nu_corrigido

    def CoeficienteDeTransferenciaDeCalor(self, densidade, viscosidade, pr_medio, pr_superficie, numero_fileiras, condutividade_termica):
        """
        Calculates the heat transfer coefficient.
        """
        _, nu_corrigido = self.NumeroDeNusselt(densidade, viscosidade, pr_medio, pr_superficie, numero_fileiras)
        coeficiente_h = nu_corrigido * (condutividade_termica / self.diametro)
        return coeficiente_h

    def VazaoMassica(self, densidade, numero_de_celulas_transversal):
        """
        Calculates mass flow rate.
        """
        return densidade * self.v_in * numero_de_celulas_transversal * self.st * self.l

    def TransferenciaDeCalor(self, densidade, numero_de_celulas_transversal, numero_de_celulas_longitudinal, temp_super, temp_entr, viscosidade, pr_medio, pr_superficie, condutividade_termica, calor_especifico):
        """
        Calculates heat transfer using various inputs and corrections.
        """
        vazao_massica = self.VazaoMassica(densidade, numero_de_celulas_transversal)
        coeficiente_h = self.CoeficienteDeTransferenciaDeCalor(densidade, viscosidade, pr_medio, pr_superficie, numero_de_celulas_longitudinal, condutividade_termica)
        
        # Total heat transfer area
        area_transferencia_calor = numero_de_celulas_longitudinal * numero_de_celulas_transversal * self.diametro * self.l

        # Logarithmic mean temperature difference
        temp_saida = temp_super - (temp_super - temp_entr) * exp(-(area_transferencia_calor * coeficiente_h) / (vazao_massica * calor_especifico))
        
        diferenca_media_log_temp = ((temp_super - temp_saida) - (temp_super - temp_entr)) / (log((temp_super - temp_saida) / (temp_super - temp_entr)))

        transf_de_calor = coeficiente_h * area_transferencia_calor * diferenca_media_log_temp

        return temp_saida, diferenca_media_log_temp, transf_de_calor

# Testing the code

densidade = 1.145
viscosidade = 1.895e-5
pr_medio = 0.7268
pr_superficie = 0.7255
numero_fileiras = 6
condutividade_termica = 0.02625
numero_de_celulas_transversal = 5
temp_super = 40
temp_entr = 30
calor_especifico = 1007

# seg = Segmento(0.018, 0.001, 0.001, 0.065, 5)
# v_max = seg.VelocidadeMaxima()
# re = seg.NumeroDeReynolds(densidade, viscosidade)
# nu, nu_corrigido = seg.NumeroDeNusselt(densidade, viscosidade, pr_medio, pr_superficie, numero_fileiras)
# coeficiente_h = seg.CoeficienteDeTransferenciaDeCalor(densidade, viscosidade, pr_medio, pr_superficie, numero_fileiras, condutividade_termica)
# vazao_massica = seg.VazaoMassica(densidade, numero_de_celulas_transversal)
# temp_saida, diferenca_media_log_temp, transf_de_calor = seg.TransferenciaDeCalor(densidade, numero_de_celulas_transversal, numero_fileiras, temp_super, temp_entr, viscosidade, pr_medio, pr_superficie, condutividade_termica, calor_especifico)

espacamento_transversal = np.arange(0.001, 0.021, 0.001)
espacamento_longitudinal = np.arange(0.001, 0.021, 0.001)
transf_de_cal = []

for espacamento_tr in espacamento_transversal:
  for espacamento_long in espacamento_longitudinal:

    seg = Segmento(0.018, espacamento_tr, espacamento_long, 0.065, 5)
    v_max = seg.VelocidadeMaxima()
    re = seg.NumeroDeReynolds(densidade, viscosidade)
    nu, nu_corrigido = seg.NumeroDeNusselt(densidade, viscosidade, pr_medio, pr_superficie, numero_fileiras)
    coeficiente_h = seg.CoeficienteDeTransferenciaDeCalor(densidade, viscosidade, pr_medio, pr_superficie, numero_fileiras, condutividade_termica)
    vazao_massica = seg.VazaoMassica(densidade, numero_de_celulas_transversal)
    temp_saida, diferenca_media_log_temp, transf_de_calor = seg.TransferenciaDeCalor(densidade, numero_de_celulas_transversal, numero_fileiras, temp_super, temp_entr, viscosidade, pr_medio, pr_superficie, condutividade_termica, calor_especifico)
    transf_de_cal.append(transf_de_calor)


data = {
    'Espacamento Transversal': np.repeat(espacamento_transversal, len(espacamento_longitudinal)),
    'Espacamento Longitudinal': np.tile(espacamento_longitudinal, len(espacamento_transversal)),
    'Transf de Cal': transf_de_cal
}

df = pd.DataFrame(data)

# Criando o gráfico de dispersão
fig = px.scatter(df, x='Espacamento Transversal', y='Espacamento Longitudinal', color='Transf de Cal',
                 color_continuous_scale='Viridis', title='Gráfico de Dispersão de Transferência de Calor')

fig.show()
