# Phy607_Project3
### Installing the Package:

To install this package, follow this procedure: 

```
cd /path/you/want/repo 
```

then clone the repo and cd into it:

```
/git clone git@github.com:jrrice5908/PHY607_Project3.git 
cd PHY607_Project3
```

Now create the conda environment:


```
conda create -n "PHY607" python=3.12
```

Enter the environment:

```
conda activate PHY607
```
then install the package:

```
pip install -e .
```

### Running the Main Script:
Finally, to run the main script, you'll need to make it an executable:

```bash
chmod +x main.sh
./main.sh
```
