from flask import Flask, render_template, request
import numpy as np
import sympy as sp

app = Flask(__name__, static_url_path='/static')


@app.route('/')
def main_menu():
    return render_template('main_menu.html')


@app.route('/continuous_stirred_reactor')
def continuous_stirred_reactor():
    return render_template('index.html')

@app.route('/calculate_csrt', methods=['POST'])
def calculate_csrt():

    flowrate = float(request.form['flowrate'])
    conversion = float(request.form['conversion'])
    rate_constant = float(request.form['rate_constant'])
    concentration = float(request.form['concentration'])
    alpha = float(request.form['alpha'])

    V = (flowrate * conversion) / ((rate_constant) * ((concentration) * (1 - conversion))**alpha)
    T = V / (flowrate/concentration)
    
    return render_template('result.html', V=V , T=T)

@app.route('/calculateBatch', methods=['POST'])
def calculate_batch():

    concentration = float(request.form['concentration'])
    conversion = float(request.form['conversion'])
    rate_constant = float(request.form['rate_constant'])
    alpha = float(request.form['alpha'])

    if(alpha==1):
        tag=np.log(1-conversion)/(-1*rate_constant)
    else:
        tag=((concentration*(1-conversion))**(1-alpha)-concentration**(1-alpha))/(-1*rate_constant*(1-alpha))
    
    return render_template('result_batch.html', Tag=tag, C=conversion)

@app.route('/pbr_liquid')
def pbr_liquid():
    return render_template('index_pbr_liquid.html')

@app.route('/calculatePbr_liquid', methods=['POST'])
def calculate_Pbr_liquid():

    concentration = float(request.form['concentration'])
    conversion = float(request.form['conversion'])
    rate_constant = float(request.form['rate_constant'])
    alpha = float(request.form['alpha'])
    flow_rate = float(request.form['flow_rate'])

    x = sp.symbols('x')
    f = 1/((1-x)**alpha)
    F_x = sp.integrate(f, (x, 0, conversion))

    cat_weight = (flow_rate * F_x) / (rate_constant*(concentration**alpha))

    Ra_dash= -1*(rate_constant*(concentration*(1-conversion))**(alpha))
    
    return render_template('result_pbr_liquid.html', Ra_dash=Ra_dash, Cat_weight=cat_weight)

@app.route('/pbr_gas')
def pbr_gas():
    return render_template('index_pbr_gas.html') 

@app.route('/calculatePbr_gas', methods=['POST'])
def calculate_Pbr_gas():

    concentration = float(request.form['concentration'])
    conversion = float(request.form['conversion'])
    rate_constant = float(request.form['rate_constant'])
    alpha = float(request.form['alpha'])
    a = float(request.form['coff_a'])
    b = float(request.form['coff_b'])
    c = float(request.form['coff_c'])
    d = float(request.form['coff_d'])
    reactor_length = float(request.form['reactor_length'])
    flow_rate = float(request.form['flow_rate'])
    gas_density = float(request.form['gas_density'])
    # superficial_velocity = float(request.form['superficial_velocity'])
    porosity = float(request.form['porosity'])
    particle_diameter = float(request.form['particle_diameter'])
    gas_viscosity = float(request.form['gas_viscosity'])
    epsilon = c+d-(a+b)
    pressure_not = float(request.form['pressure_not'])
    # volume_not = float(request.form['volume_not'])
    superficial_mass_velocity = float(request.form['superficial_mass_velocity'])
    gc = float(request.form['Gc'])

    beta_not_term1 = ((superficial_mass_velocity/(gas_density*gc))*(1-porosity)/(particle_diameter*(porosity**3)))
    beta_not_term2 =  (150*(1-porosity)*gas_viscosity/particle_diameter)+1.75*(superficial_mass_velocity)
    beta_not = beta_not_term1*beta_not_term2/(144*14.7)
    print(beta_not)

    p=(1-(2*beta_not*reactor_length*(1+epsilon))/pressure_not)**(0.5)

    ca = concentration*((1-conversion)/(1+epsilon*conversion))*p
    Ra_dash = -1*(rate_constant*(ca**alpha))
    x = sp.symbols('x')
    f = ((1+epsilon*x)/(1-x))**alpha
    F_x = sp.integrate(f, (x, 0, conversion))

    cat_weight = (flow_rate * F_x) / (rate_constant*(concentration**alpha)*(p**alpha))

    return render_template('result_pbr_gas.html', Ra_dash=Ra_dash, Cat_weight=cat_weight)


@app.route('/pfr_liquid')
def pfr_liquid():
    return render_template('index_pfr_liquid.html')

@app.route('/calculatePfr_liquid', methods=['POST'])
def calculate_Pfr_liquid():

    concentration = float(request.form['concentration'])
    conversion = float(request.form['conversion'])
    rate_constant = float(request.form['rate_constant'])
    alpha = float(request.form['alpha'])
    flow_rate = float(request.form['flow_rate'])

    Ra= -1*rate_constant*(concentration*(1-conversion))**(alpha)

    x = sp.symbols('x')
    f = 1/((1-x)**alpha)
    F_x = sp.integrate(f, (x, 0, conversion))

    Volume = (flow_rate *  F_x)/(rate_constant*(concentration**alpha))
    
    return render_template('result_pfr_liquid.html', Ra=Ra, Volume=Volume)

@app.route('/pfr_gas')
def pfr_gas():
    return render_template('index_pfr_gas.html') 

@app.route('/calculatePfr_gas', methods=['POST'])
def calculate_Pfr_gas():

    concentration = float(request.form['concentration'])
    conversion = float(request.form['conversion'])
    rate_constant = float(request.form['rate_constant'])
    alpha = float(request.form['alpha'])
    a = float(request.form['coff_a'])
    b = float(request.form['coff_b'])
    c = float(request.form['coff_c'])
    d = float(request.form['coff_d'])
    flow_rate = float(request.form['flow_rate'])
    epsilon = c+d-(a+b)
    input_pressure_not = float(request.form['input_pressure_not'])
    output_pressure_not = float(request.form['output_pressure_not'])

    ca = concentration*(1-conversion)*output_pressure_not/((1+epsilon*conversion)*input_pressure_not)
    Ra = -1*(rate_constant*(ca**alpha))

    x = sp.symbols('x')
    f = ((1+epsilon*x)/(1-x))**alpha
    F_x = sp.integrate(f, (x, 0, conversion))

    Volume = ((flow_rate * F_x)*(input_pressure_not)**alpha)/(rate_constant*(concentration*output_pressure_not)**alpha)

    return render_template('result_pfr_gas.html', Ra=Ra, Volume=Volume)


@app.route('/batch_reactor')
def batch_reactor():
    return render_template('index_batch.html')

@app.route('/plug_flow_reactor')
def plug_flow_reactor():
    return render_template('index_pfr.html')

@app.route('/packed_bed_reactor')
def packed_bed_reactor():
    return render_template('index_pbr.html')

if __name__ == '__main__':
    app.run(debug=True)
