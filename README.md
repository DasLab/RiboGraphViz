# RiboGraphViz

Visualize global properties of large RNAs using force-directed layout from GraphViz.

Authors:
* Hannah Wayment-Steele
* Credit to utils from Rhiju Das' ToyFold-1D code.
![](images/MS2_example.png)
*Above*: MS2 bacteriophage genome structure, colored by prob(unpaired), calculated in EternaFold.

*Below*: Visualizing the MFE structure and p(unpaired) of an mRNA for eGFP at increasing temperatures in Vienna.
![](images/melting_eGFP_mRNA.gif)

Note: Not intended for detailed layouts -- loops may switch orientation in z-axis.

To set up:
```
sudo pip install -r requirements.txt
sudo python setup.py install
```

You'll need to use Python3; Python2 won't work. 

*Tips for Mac users*: if you're stuck with Python 2, setting up a virtual environment (e.g., with `conda environment`) to install Python3 might be a good choice. You may also need to install `pygraphviz` -- which does not work with conda, but does work with `brew` and `pip install`, as noted [here](https://pygraphviz.github.io/documentation/stable/install.html). 

See `examples.ipynb` for example usage.

![](images/RGV_example_colorings.png)

*Below*: Vienna MFE structures of 20 randomly-generated RNAs.
![](images/multiple_struct_example.png)
