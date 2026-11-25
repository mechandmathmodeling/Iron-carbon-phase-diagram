import os
import sys
import time
from playwright.sync_api import sync_playwright

APP_URL = os.environ.get("STREAMLIT_APP_URL")
if not APP_URL:
    print("Error: STREAMLIT_APP_URL environment variable is not set.")
    sys.exit(1)

with sync_playwright() as p:
    browser = p.chromium.launch(headless=True)
    page = browser.new_page()
    print(f"Navigating to {APP_URL}...")
    page.goto(APP_URL, wait_until="networkidle", timeout=60000)

    # Check if the app is asleep and click the wake-up button if present
    wake_button = page.locator("button:has-text('Yes, get this app back up!')")
    if wake_button.is_visible():
        print("App was asleep. Clicking wake-up button...")
        wake_button.click()
        # Wait 45 seconds for the container to boot and connect
        time.sleep(45)
    else:
        print("App is already awake.")
        # Keep the WebSocket connection open briefly so Streamlit logs the visit
        time.sleep(10)

    browser.close()
    print("Done.")
