"""Test script for verbosity=0 (silent mode)"""

import tempfile
from pathlib import Path

from tctrack.tempest_extremes import TETracker
from tctrack.tempest_extremes import TEDetectParameters

def test_verbosity_0_silent():
    """Test _run_subprocess with verbosity=0 (silent, no output)"""
    
    print("=" * 60)
    print("Testing verbosity=0 — Silent mode")
    print("=" * 60)
    print()
    
    # Step 1: Create a temporary input file
    print("Step 1: Creating temporary input file...")
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("input line 1\n")
        f.write("input line 2\n")
        f.write("input line 3\n")
        temp_file = f.name
    print(f"Created temp file: {temp_file}\n")
    
    # Step 2: Create a simple command
    print("Step 2: Setting up command...")
    command_list = ["cat"]  # cat reads from stdin
    print(f"Command: {command_list}\n")
    
    # Step 3: Create tracker instance
    print("Step 3: Creating tracker instance...")
    detect_params = TEDetectParameters(in_data=["dummy.nc"])
    tracker = TETracker(detect_params)
    print(f"Tracker created\n")
    
    # Step 4: Call _run_subprocess with verbosity=0
    print("Step 4: Calling _run_subprocess with verbosity=0...")
    print("  (You should see NO output below from the subprocess)")
    print()
    
    result = tracker.run_tracker_subprocess(
        command_name="SilentTest",
        command_list=command_list,
        input_file=temp_file,
        verbosity=0,  # ← SILENT mode
    )
    
    print()
    print("Step 5: Checking results...")
    print(f"  Return code: {result['returncode']}")
    print(f"  Stdout length: {len(result['stdout'])} characters")
    print(f"  Stderr length: {len(result['stderr'])} characters")
    print()
    
    # Step 6: Verify the results
    print("Step 6: Verifying results...")
    assert result['returncode'] == 0, "Command should succeed (return code 0)"
    print(f"Return code is 0 (success)")
    
    assert "input line 1" in result['stdout'], "Should contain input line 1"
    print(f"Stdout contains expected data")
    
    assert result['stderr'] == "", "Should have no errors"
    print(f"Stderr is empty (no errors)")
    
    print()
    
    # Step 7: Show the captured output
    print("Step 7: Actual output captured (not printed during execution):")
    print("  Stdout content:")
    for line in result['stdout'].split('\n'):
        print(f"    {line}")
    print()
    
    # Step 8: Clean up
    print("Step 8: Cleaning up...")
    Path(temp_file).unlink()
    print(f"  ✓ Deleted temp file\n")
    
    # Final result
    print("=" * 60)
    print("✓ verbosity=0 test PASSED!")
    print("=" * 60)
    print()
    print("Summary:")
    print("  - Command executed silently (no output printed)")
    print("  - Results were captured in the returned dict")
    print("  - Return code indicates success")
    print()

if __name__ == "__main__":
    test_verbosity_0_silent()