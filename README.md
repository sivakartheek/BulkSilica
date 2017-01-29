# Bulk silica model
## An interatomic potentials based Metropolis Monte Carlo model for investigating the structural, thermodynamic and bulk properties of silicon dioxide (silica).

Using this program (main_silica.cpp), simulations were performed with 64-SiO2 groups (64 silicon and 128 oxygen atoms) at different temperatures and pressures. All the simulations start from the crystal structure of SiO2 (beta-cristobalite) to avoid atoms being close to each other. Each simulation consisted of an initial equilibration of 1.5 million moves followed by additional 0.5 million production moves. The red line represents the production region, where averages computed.

![alt tag](https://raw.githubusercontent.com/NaveenKaliannan/BulkSilica/master/output/graph/Density_.png)

![alt tag](https://raw.githubusercontent.com/NaveenKaliannan/BulkSilica/master/output/graph/Energy.png)

## Density
The first property we discuss here is the density of silica glass, which depends on the temperature and pressure. The density of silica glass obtained from the simulations under different temperatures and pressures is shown in the below figure. Clearly higher the
temperature, lower the density of SiO2. We have tted straight lines to describe the temperature-density behaviour for different pressures. It is likely that the singular behaviour at 10 GPa is due to insuffcient equilibration.

![alt tag](https://raw.githubusercontent.com/NaveenKaliannan/BulkSilica/master/output/graph/density.png)

## Structural properties
In order to study the silica glass structure, and its variation with temperature and pressure, we calculate the total and partial RDFs under different temperatures and pressures as shown in the below Figures. The first and second peak positions of Si-Si, O-O and Si-O RDFs correspond to the Si-Si, O-O and Si-O mean rst and second nearest neighbour distances, respectively. In the total RDF, the frst, second and third peak positions correspond to the Si-O, O-O and Si-Si mean frst nearest neighbour distances, respectively.

![alt tag](https://raw.githubusercontent.com/NaveenKaliannan/BulkSilica/master/output/graph/total.png)

![alt tag](https://raw.githubusercontent.com/NaveenKaliannan/BulkSilica/master/output/graph/si.png)

![alt tag](https://raw.githubusercontent.com/NaveenKaliannan/BulkSilica/master/output/graph/sio.png)

![alt tag](https://raw.githubusercontent.com/NaveenKaliannan/BulkSilica/master/output/graph/o.png)




