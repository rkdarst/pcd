
docs-clean:
	rm -r doc/api/ || true
	rm -r doc/_build/ || true

docs:
	sphinx-apidoc -P -e -o doc/api/ -H "``pcd`` API" ../pcd/ tests lazy.py
#	cd .. && sphinx-apidoc -P -e -o pcd/doc/api/ pcd/ pcd/tests pcd/lazy.py
	cd doc && PYTHONPATH=..:${PYTHONPATH} make html
#	cd doc && PYTHONPATH=..:$PYTHONPATH make html
