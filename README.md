# RiboGraphViz

Authors:
* Hannah Wayment-Steele
* Credit to utils from Rhiju Das' ToyFold-1D code.

Visualize global properties of large RNAs, using force-directed layout from GraphViz.

![](images/MS2_example.png)

*Above*: MS2 bacteriophage genome structure, colored by prob(unpaired), calculated in EternaFold.

*Below*: Visualizing the MFE structure and p(unpaired) of an mRNA for eGFP, at increasing temperatures in Vienna.
![](images/melting_eGFP_mRNA.gif)

Note: Not intended for detailed layouts -- loops may switch orientation in z-axis.

To set up:
```
pip install -r requirements.txt
python setup.py install
```

See `examples.ipynb` for example usage.

![](images/RGV_example_colorings.png)

*Below*: Vienna MFE structures of 20 randomly-generated RNAs.
![](images/multiple_struct_example.png)
