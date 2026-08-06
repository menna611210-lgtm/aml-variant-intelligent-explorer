import streamlit as st

from components.styles import load_css
from pages.home import show

st.set_page_config(
    page_title="Genova",
    page_icon="🧬",
    layout="wide"
)

load_css()

show()