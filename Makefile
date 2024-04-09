.PHONY: push jupyter

ifdef m
push:
	git add -A
	git commit -am "$(m)"
	git push
else
push:
	$(error Please provide a commit message, e.g., make push m="your message here")
endif

jupyter:
	jupyter notebook --ip 0.0.0.0 --port 8888 --no-browser --allow-root

jupyterlocal:
	sudo jupyter notebook --allow-root

