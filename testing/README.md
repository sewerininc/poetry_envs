To create a env:
1. mkdir project name
2. poetry init
3. poetry add ipykernel
4. poetry install
   # This next line will create the kernel to be used in jupyter notebook if you update the poetry env you will have to make a this line again (the same name will overwrite)
6. poetry python -m ipykernel install --user --name=my_poetry_env
7. jupyter notebook
8. select your env for jupyter to be used - this is changed on the kernel button
