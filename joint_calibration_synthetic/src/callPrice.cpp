#pragma once

#include <ql/quantlib.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace qe {

double hestonCallPrice(
    double spot,
    double strike,
    double maturityYears,
    double r,
    double q,
    double v0,
    double kappa,
    double theta,
    double xi,
    double rho,
    const QuantLib::Date& evalDate) {

    using namespace QuantLib;

    if (!(spot > 0.0) || !(strike > 0.0)) {
        throw std::invalid_argument("spot and strike must be > 0");
    }
    if (!(maturityYears > 0.0)) {
        throw std::invalid_argument("maturityYears must be > 0");
    }
    if (!(v0 >= 0.0) || !(kappa > 0.0) || !(theta >= 0.0) || !(xi >= 0.0)) {
        throw std::invalid_argument("invalid Heston variance parameters");
    }

    const double clampedRho = std::clamp(rho, -0.999, 0.999);

    Date today = evalDate;
    if (today == Date()) {
        today = Settings::instance().evaluationDate();
        if (today == Date())
            today = Date::todaysDate();
    }
    if (today == Date())
        throw std::runtime_error("invalid evaluation date");
    Settings::instance().evaluationDate() = today;

    DayCounter dc = Actual365Fixed();
    Calendar cal = TARGET();

    const Integer nDays = std::max<Integer>(1, static_cast<Integer>(std::lround(365.0 * maturityYears)));
    const Date maturityDate = cal.advance(today, nDays, Days);

    Handle<Quote> spotHandle(ext::make_shared<SimpleQuote>(spot));
    Handle<YieldTermStructure> rTS(ext::make_shared<FlatForward>(today, r, dc));
    Handle<YieldTermStructure> qTS(ext::make_shared<FlatForward>(today, q, dc));

    auto process = ext::make_shared<HestonProcess>(
        rTS, qTS, spotHandle,
        v0, kappa, theta, xi, clampedRho);
    auto model = ext::make_shared<HestonModel>(process);
    auto engine = ext::make_shared<AnalyticHestonEngine>(model);

    auto payoff = ext::make_shared<PlainVanillaPayoff>(Option::Call, strike);
    auto exercise = ext::make_shared<EuropeanExercise>(maturityDate);
    VanillaOption option(payoff, exercise);
    option.setPricingEngine(engine);

    return option.NPV();
}

}
