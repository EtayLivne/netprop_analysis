from time import sleep
from typing import Union

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager

strings = Union[str, list[str]]


def _listify(elm) -> list:
    if isinstance(elm, list):
        return elm
    return [elm]


def _build_query(path_to_genes_list: str):
    _webgestalt_url = "http://www.webgestalt.org/option.php?organism=hsapiens&enrich_method=ORA&fdr_m\
    ethod=BY&enriched_database_category=geneontology&enriched_database_name=Bio\
    logical_Process_noRedundant&sig_method=top&sig_value=10&fdr_method=BH&max_num=200&id_ty\
    pe=entrezgene&gene_list={}&ref_set=genome"
    with open(path_to_genes_list, "r") as handler:
        content = handler.read()
    return _webgestalt_url.format(content).replace("\n", '0%A')


def _submit_query(driver: webdriver.WebDriver) -> None:
    try:
        submit_button = WebDriverWait(driver, 20).until(EC.presence_of_element_located((By.XPATH,
                                                                                        '//*[@id="wgForm"]/button')))
    except:
        print("ERROR!  **** webgestalt took too long to load ****")
        driver.quit()
        exit()

    sleep(2)
    submit_button.click()

    sleep(3)
    driver.switch_to.window(driver.window_handles[1])


def _download_query_result(driver: webdriver.Webdriver) -> None:
    try:
        download = WebDriverWait(driver, 30).until(
            EC.presence_of_element_located((By.XPATH, '//*[@id="app"]/div[2]/div/div[2]/div/a'))
        )
    except:
        print("webstalt took too long to process")
        driver.quit()
        exit()

    download.click()


def query_webgestalt(paths_to_genes_lists: strings, download_path_root: str) -> None:

    for gene_list_path in _listify(paths_to_genes_lists):

        query_url = _build_query(paths_to_genes_lists)
        chrome_options = webdriver.ChromeOptions()
        prefs = {'download.default_directory': download_path}
        chrome_options.add_experimental_option('prefs', prefs)
        driver = webdriver.Chrome(ChromeDriverManager().install(), chrome_options=chrome_options)
        driver.maximize_window()
        driver.get(query_url)


# //*[contains(text(), 'Result Download')]