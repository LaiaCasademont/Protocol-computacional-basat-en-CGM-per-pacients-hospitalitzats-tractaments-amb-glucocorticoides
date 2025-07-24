# Protocol-computacional-basat-en-CGM-per-pacients-hospitalitzats-tractaments-amb-glucocorticoides
Aquest Treball de Fi de Grau proposa un protocol computacional per optimitzar el control glucèmic en pacients hospitalitzats no crítics tractats amb glucocorticoides, mitjançant el monitoratge continu de glucosa (CGM). A partir d’un simulador basat en el model de Hovorka, desenvolupat per investigadors de la Càtedra UdG–Dexcom, s’han aplicat modificacions específiques per estudiar l’efecte dels glucocorticoides i validar un algoritme terapèutic amb alarmes predictives.


## Estructura del repositori
README.md: aquest fitxer, amb la descripció del projecte i instruccions d’ús.

**->**code/: carpeta que conté els fitxers .m creats o modificats durant el projecte:
_Aquest TFG s’ha desenvolupat a través de nombroses versions internes. En aquest repositori es publica únicament la versió final.
El procés complet, les iteracions prèvies i l’evolució del protocol es troben detallats a la memòria del TFG._
_Si considereu d’utilitat accedir a versions anteriors o a altres etapes del desenvolupament, no dubteu a contactar-me per correu electrònic:
laiacasademontt@gmail.com_

**-** main.m: codi principal d'execució de la simulació (on s'han implentat les modifiacions per simular la hiperglucèmia produïda pels glucocorticoides.

**-** run_mycontroller.m: implementació del protocol terapèutic final desenvolupat.

**-** init_mycontroller.m: configuració inicial del controlador amb les variables inicialitzades.

**-** generate_glucose.m: funció per simular glucosa capil·lar a partir del CGM.


**->**outputs/: resultats rellevants generats per la simulació.
exemple_resultats.mat: dades de sortida simulades.

**->**docs/: documents de suport.
A_simulator_with_realistic_and_challenging_scenarios_for_virtual_T1D_patients_undergoing_CSII_and_MDI_therapy.pdf: paper on es desenvolupa el simulador original al qual s'han aplicat els canvis.
Memòria_TFG_Laia_Casademont.pdf: versió final de la memòria del TFG.

## Requeriments d'instal·lació 
- **MATLAB R2024a** o superior  
- Toolboxes recomanats:
  - Signal Processing
  - Statistics and Machine Learning
    
**AVÍS DE PROPIETAT INTEL·LECTUAL!**

**El simulador base emprat en aquest projecte és propietat de la Càtedra UdG–Dexcom. fou desenvolupat en el marc de la Tesi Doctoral de l'Ernesto Estremera, investigador de la Càtedra i professor de la UdG. No s’inclou en aquest repositori per motius de drets d’autor.
Per accedir al simulador original, cal sol·licitar permís directament a través de la Càtedra UdG–Dexcom**
