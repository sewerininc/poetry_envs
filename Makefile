.PHONY: push

ifeq ($(MESSAGE),)
$(error Please provide a commit message, e.g., make push MESSAGE="your message here")
endif

push:
	git add -A
	git commit -am "$(MESSAGE)"
	git push

