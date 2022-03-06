# Growth curve fittR: *Easy growth curve curve fitting for experimental biologists*

Scott H. Saunders, Distinguished Fellow UTSW, ![lab website]()

---

* ![link to app]()
* ![link to preprint / pdf]()

## Motivation

One of the most common experiments in microbiology is the growth curve. Typically, absorbance measurements (wavelength = 600 nm) are taken over time, which correspond to cell density in the liquid medium. As cells divide to reproduce, their overall growth reflects an exponential process, until they run out of food and stop growing. These growth curves are often used as qualitative data to show that cells grew more slowly or to a lower density in one condition compared to another. However, there are many ways that growth curves can be understood quantitatively. Microbiologists have long appreciated this fact, and generations of trainees have calculated parameters by finding the slope of a growth curve on a log scale. This project aims to modernize this tradition, by empowering users to rapidly quantify parameters and uncertainty from their data without having to code.

![Screenshot images of what app does]()

## Details

<p align="center">
<img src="gomp_diagram.png" width=50% height=50%>
</p>


<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=y=Aexp{\left(-exp{\left(\frac{\mu_A e}{A}(\lambda-t)+1\right)}\right)} +C" width=30% height=30%>
</p>

