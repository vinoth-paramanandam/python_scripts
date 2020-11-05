# from selenium import webdriver

# driver = webdriver.Chrome()

# driver.get("https://netaccess.iitm.ac.in")
# assert "iitm" in driver.title

# username = driver.find_element_by_id("username")
# username.clear()
# username.send_keys("ae15d018")

# password = driver.find_element_by_id("password")
# password.clear()
# password.send_keys("7mGJyPmJ")

# driver.find_element_by_id("submit").click()
from selenium import webdriver

webdriver.Chrome().get("https://netaccess.iitm.ac.in")