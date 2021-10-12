

clean-install:
	pip uninstall flares_utility
	pip install . -r requirements.txt

install:
	pip install . -r requirements.txt
