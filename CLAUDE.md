# Test design

1. When running a specific test, use `pixi run run-a-test test_name`, where `test_name` is the name of the test function to run (it is specified as the argument for pytest's `-k`). This allows for quick iteration on a specific test without running the full test suite.
2. To test for changes more broadly (but still relatively quickly), use `pixi run run-quick-tests`
3. To do the final expensive tests, use `pixi run run-tests-at-cmr`
4. Add tests to the appropriate test files in the `test` directory. For example for changes to the coassemble subcommand, add tests to `test/test_coassemble.py`. For changes to the `cluster_graph.py` script, add tests to `test/test_cluster_graph.py`. For changes that run expensive operations (e.g. full-sample assembly), add tests to `test/test_manual.py`.
5. Update the documentation in docs/ to include any relevant information about the new tool, including usage instructions and examples.
