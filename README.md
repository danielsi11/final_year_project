# **Final Year Project**

## **Section_9.1**
This folder contains all the source codes used to simulate a single gateway system with a single group cluster. Inside the folder, the following subfolders contains the necessary files used to produce the respective figures in the report:

### **Figure_6** 
In this subfolder, the topology of a single beam cluster with 7 beams will be produced by running the "topology.m" file. By changing the formation of the beam centres, a different pattern can be obtained. By changing the height or the 3 dB angle, a different radius for the single beam coverage can be obtained.

### **Figure_7**
In this subfolder, the system under both perfect and imperfect CSIT will be implemented and the results will be plotted on the same figure by running the "main.m" file. There are also MATLAB functions created to implement the optimisation problem for the 1-layer RS and NoRS transmission schemes.

***

## **Section_9.2**
The folder contains all the source codes used to simulate a single gateway system with multiple group clusters. Inside the folder, the following subfolders contains the necessary files used to produce the respective figures in the report:

### **Figure_8**
In this subfolder, the topology of a single beam cluster with 9 beams and multiple group clusters with 3 groups each will be produced by running the "topology.m" file. By changing the number of beam boundaries to be plotted with the same colour, a different group cluster size can be obtained.

### **Figure_9**
In this subfolder, the system under perfect CSIT will be implemented and the results will be plotted on the same figure by running the "main.m" file. In addition to the 1-layer RS and NoRS, there are also MATLAB functions created to implement the optimisation problem for the RS per cluster, HRS and GRS transmission schemes.

### **Figure_10**
In this subfolder, the same system but under imperfect CSIT will be implemented and the results will be plotted on the same figure by running the "main.m" file. The MATLAB functions created to implement the 1-layer RS, NoRS, RS per cluster, HRS and GRS transmission schemes were adjusted to carry out the SAA method to form the deterministic version of the optimisation problem.

### **Figure_11**
In this subfolder, the topology of a single beam cluster with 9 beams and multiple group clusters with 3 groups each but with a different group clustering will be produced by running the "topology.m" file. By changing the sequence of the beam centres, a different group clustering can be obtained.

###  **Figure_12**
In this subfolder, the same system but with a different group clustering under both perfect and imperfect CSIT will be implemented and the results will be plotted on the same figure by running the "main.m" file. The MATLAB functions created to implement the 1-layer RS, NoRS, RS per cluster and HRS transmission schemes were adjusted to include the change in the group clustering.

***

## **Section_9.3**
The folder contains all the source codes used to simulate a multiple gateway system. Inside the folder, the following subfolders contains the necessary files used to produce the respective figures in the report:

### **Figure_13**
In this subfolder, the topology of multiple beam clusters with 3 beams each will be produced by running the "topology.m" file. By changing the sequence of the beam centres and the number of beam boundaries to be plotted with the same colour, a different cluster of groups served by the same gateway can be obtained.

### **Figure_14**
In this subfolder, the system without feeder link interference and noise under both perfect and imperfect CSIT will be implemented and the results will be plotted on the same figure by running the "main.m" file. The MATLAB functions created to implement the optimisation problem for the 1-layer RS and NoRS transmission schemes were adjusted to be implemented at each gateway.

### **Figure_15**
In this subfolder, the same system but with feeder link interference and noise under both perfect and imperfect CSIT will be implemented and the results will be plotted on the same figure by running the "main.m" file. The MATLAB functions created to implement the 1-layer RS and NoRS transmission schemes at each gateway were adjusted to include the feeder link interference and noise.

### **Figure_16**
In this subfolder, the same system but with on-board processing together with feeder link interference and noise under both perfect and imperfect CSIT will be implemented and the results will be plotted on the same figure by running the "main.m" file. The MATLAB functions created to implement the 1-layer RS and NoRS transmission schemes at each gateway were adjusted to include the on-board processing and the two-stage precoding.