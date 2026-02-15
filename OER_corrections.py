import scipy.constants as sc
import numpy as np
from typing import (List,)
import tomllib

boltzmann_constant = sc.physical_constants["Boltzmann constant in eV/K"][0] # eV K-1
planck_constant = sc.physical_constants["Planck constant in eV/Hz"][0] # eV s
speed_of_light = sc.speed_of_light/sc.centi # the speed of light in cm/s
gas_constant_R = boltzmann_constant # eV K-1

# NOTE: gas_constant_R assumes this unit because it is divided by Avogadro, ie the resulting energy is per molecule.

class Corrections:

    def __init__(self, adsorbate_name: str, frequencies: List, temperature: float):

        if not isinstance(adsorbate_name, str):
            raise TypeError("The name of the adsorbed specimen must be a string.")
        elif len(adsorbate_name) == 0:
            raise ValueError("The name of the adsorbate must contain at least one character.")
        if not isinstance(frequencies, List):
            raise TypeError("Frequencies must be of type List.")
        elif len(frequencies) == 0:
            raise ValueError("The list of frequencies must contain at least one value.")
        for freq in frequencies:
            if not isinstance(freq, float):
                raise TypeError(f"Frequencies must be floating point values. Invalid element: {freq}.")
        if not isinstance(temperature, float):
            raise TypeError("Temperature must be a floating point value.")
        elif temperature <= 0:
            raise ValueError("Temperature must be a positive value.")
        
        self.name = adsorbate_name
        self.frequencies = frequencies
        self.temperature = temperature
        self.ZPE = self.calculate_ZPE(self.frequencies)
        self.S_vib = self.calculate_Svib(self.frequencies, self.temperature)
        self.deltaU = self.calculate_deltaU(self.frequencies, self.temperature)
        print(f"Specimen: {self.name}")
        print(f"ZPE: {self.ZPE} eV\t\tS_vib: {self.S_vib} eV\t\tT*S_vib: {self.temperature*self.S_vib}\t\tdeltaU 0->T: {self.deltaU} eV")
        print(f"Total correction: {self.ZPE - self.temperature*self.S_vib + self.deltaU} eV")

    @staticmethod
    def calculate_ZPE(frequency_list: List) -> float:
        ''' 
            Evaluates the zero-point energy (ZPE, in eV units) from a list of vibrational frequencies (in cm-1).

            Parameters
            ----------
            frequency_list: 
                List of vibrational frequencies (in cm-1).

            Returns
            -------
            float
                ZPE correction in eV.
        '''

        return 0.5*np.sum([planck_constant*speed_of_light*freq for freq in frequency_list])

    @staticmethod
    def calculate_Svib(frequency_list: List, temperature: float):
        ''' 
            Evaluates the vibrational contribution to entropy (Svib, in eV units) from a list of vibrational frequencies (in cm-1).

            Parameters
            ----------
            frequency_list: 
                List of vibrational frequencies (in cm-1).

            Returns
            -------
            float
                Svib contribution in eV.
        '''
        hnu_over_kt = lambda f: planck_constant*speed_of_light*f/(boltzmann_constant*temperature)
        return gas_constant_R*np.sum([((planck_constant*speed_of_light*frequency)/(boltzmann_constant*temperature*(np.exp(hnu_over_kt(frequency)) - 1))) - np.log(1 - np.exp(-hnu_over_kt(frequency))) for frequency in frequency_list])

    @staticmethod
    def calculate_deltaU(frequency_list: List, temperature: float):
        ''' 
            Evaluates the internal energy variation (deltaU 0->T, in eV units) from a list of vibrational frequencies (in cm-1).

            Parameters
            ----------
            frequency_list: 
                List of vibrational frequencies (in cm-1).

            Returns
            -------
            float
                deltaU 0->T energy variation in eV.
        '''
        hnu_over_kt = lambda f: planck_constant*speed_of_light*f/(boltzmann_constant*temperature)
        return gas_constant_R*np.sum([(planck_constant*speed_of_light*frequency)/(boltzmann_constant*(np.exp(hnu_over_kt(frequency)) - 1)) for frequency in frequency_list])

if __name__ == "__main__":
    try:
        with open("data.toml", "rb") as fp:
            data = tomllib.load(fp)
            print(data)
            
    except FileNotFoundError:
        print("The .toml configuration file was not found.")
        raise
    else:
        temperature = data["temperature"]
        for i in data["specimen"]:
            Corrections(adsorbate_name=i["name"], frequencies=[j for j in i["frequencies"]], temperature=temperature)