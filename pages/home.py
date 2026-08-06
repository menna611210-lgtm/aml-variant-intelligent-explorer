import streamlit as st

def show():
    hero_section()
    action_section()
    features_section()
    footer_section()


def hero_section():
    st.markdown('''
    <div class = "hero-card">
    <h1>Genova</h1>
    <h2> AI-Powered Genomic Platform </h2>
    <p> 
       Analyze genomic datasets with a simple,
        modern and interactive interface.
    </P>
    </div>   
    ''', unsafe_allow_html=True)

def action_section():
    st.markdown('''
    <div class = "action-container">
    
      <div class = "action-card">
         <h3>Dataset Analyzer</h3>
         <p>parse and validate genomic datasets</p>
      </div>
    
      <div class = "action-card">
          <h3>VCF Analyzer</h3>
          <p>Parse and visualize genomic VCF files.</p>
      </div>
      
      <div class = "action-card">
           <h3>Demo</h3>
           <p>Explore Genova using a sample dataset.</p>
      </div>
    </div>
    ''', unsafe_allow_html=True  )
    pass

def features_section():
    pass





def footer_section():
    pass